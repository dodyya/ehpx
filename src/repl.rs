use std::io::{self, BufRead, Write};
use ehpx::Element;

// ── Lexer ─────────────────────────────────────────────────────────────────────

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Number(usize),
    Plus,
    Star,
    LBracket,
    RBracket,
    LParen,
    RParen,
    Comma,
    D,
    Eof,
}

struct Lexer {
    chars: Vec<char>,
    pos: usize,
}

impl Lexer {
    fn new(input: &str) -> Self {
        Self { chars: input.chars().collect(), pos: 0 }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn advance(&mut self) -> Option<char> {
        let c = self.chars.get(self.pos).copied();
        self.pos += 1;
        c
    }

    fn tokenize(mut self) -> Result<Vec<Token>, String> {
        let mut tokens = Vec::new();
        loop {
            while self.peek().map_or(false, |c| c.is_whitespace()) {
                self.advance();
            }
            let tok = match self.peek() {
                None => Token::Eof,
                Some('+') => { self.advance(); Token::Plus }
                Some('*') => { self.advance(); Token::Star }
                Some('[') => { self.advance(); Token::LBracket }
                Some(']') => { self.advance(); Token::RBracket }
                Some('(') => { self.advance(); Token::LParen }
                Some(')') => { self.advance(); Token::RParen }
                Some(',') => { self.advance(); Token::Comma }
                Some('d') => { self.advance(); Token::D }
                Some(c) if c.is_ascii_digit() => {
                    let mut s = String::new();
                    while self.peek().map_or(false, |ch| ch.is_ascii_digit()) {
                        s.push(self.advance().unwrap());
                    }
                    Token::Number(
                        s.parse().map_err(|_| format!("number too large: {}", s))?,
                    )
                }
                Some(c) => return Err(format!("unexpected character '{}'", c)),
            };
            let done = tok == Token::Eof;
            tokens.push(tok);
            if done { break; }
        }
        Ok(tokens)
    }
}

// ── Parser ────────────────────────────────────────────────────────────────────

struct Parser {
    tokens: Vec<Token>,
    pos: usize,
}

impl Parser {
    fn new(tokens: Vec<Token>) -> Self {
        Self { tokens, pos: 0 }
    }

    fn peek(&self) -> &Token {
        &self.tokens[self.pos]
    }

    fn advance(&mut self) -> Token {
        let tok = self.tokens[self.pos].clone();
        if self.pos + 1 < self.tokens.len() {
            self.pos += 1;
        }
        tok
    }

    fn expect(&mut self, expected: Token) -> Result<(), String> {
        let tok = self.advance();
        if tok == expected {
            Ok(())
        } else {
            Err(format!("expected {:?}, got {:?}", expected, tok))
        }
    }

    // expr = term ('+' term)*
    fn parse_expr(&mut self) -> Result<Element, String> {
        let mut acc = self.parse_term()?;
        while *self.peek() == Token::Plus {
            self.advance();
            acc = acc + self.parse_term()?;
        }
        Ok(acc)
    }

    // term = factor ('*' factor)*
    fn parse_term(&mut self) -> Result<Element, String> {
        let mut acc = self.parse_factor()?;
        while *self.peek() == Token::Star {
            self.advance();
            acc = acc * self.parse_factor()?;
        }
        Ok(acc)
    }

    // factor = '[' numlist ']' | '(' expr ')' | 'd' '(' expr ')'
    fn parse_factor(&mut self) -> Result<Element, String> {
        match self.peek().clone() {
            Token::D => {
                self.advance();
                self.expect(Token::LParen)?;
                let inner = self.parse_expr()?;
                self.expect(Token::RParen)?;
                Ok(inner.diff())
            }
            Token::LBracket => {
                self.advance();
                let mut nums: Vec<usize> = Vec::new();
                if *self.peek() != Token::RBracket {
                    match self.advance() {
                        Token::Number(n) => nums.push(n),
                        tok => return Err(format!("expected number, got {:?}", tok)),
                    }
                    while *self.peek() == Token::Comma {
                        self.advance();
                        match self.advance() {
                            Token::Number(n) => nums.push(n),
                            tok => return Err(format!("expected number after ',', got {:?}", tok)),
                        }
                    }
                }
                self.expect(Token::RBracket)?;
                Element::new(&nums)
                    .ok_or_else(|| format!("not an admissible sequence: {:?}", nums))
            }
            Token::LParen => {
                self.advance();
                let inner = self.parse_expr()?;
                self.expect(Token::RParen)?;
                Ok(inner)
            }
            tok => Err(format!("expected '[' or '(', got {:?}", tok)),
        }
    }
}

// ── Eval & REPL ───────────────────────────────────────────────────────────────

fn eval(input: &str) -> Result<Element, String> {
    let tokens = Lexer::new(input).tokenize()?;
    let mut parser = Parser::new(tokens);
    let elem = parser.parse_expr()?;
    if *parser.peek() != Token::Eof {
        return Err(format!("unexpected token: {:?}", parser.peek()));
    }
    Ok(elem)
}

pub fn run() {
    let stdin = io::stdin();
    print!("> ");
    io::stdout().flush().unwrap();
    for line in stdin.lock().lines() {
        let line = match line {
            Ok(l) => l,
            Err(_) => break,
        };
        let input = line.trim();
        if input.is_empty() {
            print!("> ");
            io::stdout().flush().unwrap();
            continue;
        }
        if input == "quit" || input == "exit" {
            break;
        }
        match eval(input) {
            Ok(elem) => println!("{}", elem),
            Err(e) => eprintln!("error: {}", e),
        }
        print!("> ");
        io::stdout().flush().unwrap();
    }
}

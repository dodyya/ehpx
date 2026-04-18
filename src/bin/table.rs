use ehpx::curtis::{CurtisTable, RenderStyle};
use std::io::{IsTerminal, Write};
use std::path::PathBuf;
use std::process::ExitCode;

fn usage() -> &'static str {
    "\
usage: table [MAX_STEM] [options]

positional:
  MAX_STEM        run Curtis algorithm through this stem (default 12)

options:
  --plain         no ANSI color; Unicode (λ, →).  Default when writing to file
                  or when stdout is not a terminal.
  --ascii         no ANSI; strict ASCII (\"l\" for lambda, \"->\").  Best for
                  portability across viewers.
  --color         force ANSI color (default when stdout is a terminal).
  --json          emit machine-readable JSON instead of the human report.
  -o, --output PATH
                  write output to PATH instead of stdout.  Turns off ANSI
                  unless --color is also given.
  -h, --help      show this help.

examples:
  table 12                                 # colored report to terminal
  table 12 --plain -o report.txt           # save .txt report
  table 12 --ascii -o report.txt           # strict ASCII .txt report
  table 12 --json -o report.json           # JSON for Python visualizer
"
}

enum Format {
    Human,
    Json,
}

struct Args {
    max_stem: usize,
    format: Format,
    style: Option<RenderStyle>, // None = auto-detect
    output: Option<PathBuf>,
}

fn parse_args() -> Result<Args, String> {
    let mut max_stem: Option<usize> = None;
    let mut format = Format::Human;
    let mut style: Option<RenderStyle> = None;
    let mut output: Option<PathBuf> = None;

    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut i = 0;
    while i < argv.len() {
        let a = &argv[i];
        match a.as_str() {
            "-h" | "--help" => {
                println!("{}", usage());
                std::process::exit(0);
            }
            "--plain" => style = Some(RenderStyle::Plain),
            "--ascii" => style = Some(RenderStyle::Ascii),
            "--color" => style = Some(RenderStyle::Ansi),
            "--json" => format = Format::Json,
            "-o" | "--output" => {
                i += 1;
                if i >= argv.len() {
                    return Err(format!("{} needs an argument", a));
                }
                output = Some(PathBuf::from(&argv[i]));
            }
            s if s.starts_with("--output=") => {
                output = Some(PathBuf::from(&s["--output=".len()..]));
            }
            s if s.starts_with("-o") && s.len() > 2 => {
                output = Some(PathBuf::from(&s[2..]));
            }
            s if !s.starts_with('-') => {
                if max_stem.is_some() {
                    return Err(format!("unexpected positional arg: {:?}", s));
                }
                max_stem = Some(
                    s.parse()
                        .map_err(|_| format!("bad MAX_STEM: {:?}", s))?,
                );
            }
            other => return Err(format!("unknown option: {:?}", other)),
        }
        i += 1;
    }

    Ok(Args {
        max_stem: max_stem.unwrap_or(12),
        format,
        style,
        output,
    })
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{}\n\n{}", e, usage());
            return ExitCode::from(2);
        }
    };

    // Auto-detect style if the user didn't force one:
    //   • writing to a file → Plain
    //   • piped stdout      → Plain
    //   • terminal stdout   → Ansi
    let style = args.style.unwrap_or_else(|| {
        if args.output.is_some() || !std::io::stdout().is_terminal() {
            RenderStyle::Plain
        } else {
            RenderStyle::Ansi
        }
    });

    eprintln!("Computing Curtis table through stem {}...", args.max_stem);
    let table = CurtisTable::compute(args.max_stem);

    let body = match args.format {
        Format::Human => table.display(args.max_stem, style),
        Format::Json => table.emit_json(args.max_stem),
    };

    match &args.output {
        Some(path) => {
            if let Err(e) = std::fs::write(path, &body) {
                eprintln!("write {}: {}", path.display(), e);
                return ExitCode::from(1);
            }
            eprintln!("wrote {}", path.display());
        }
        None => {
            let stdout = std::io::stdout();
            let mut lock = stdout.lock();
            if let Err(e) = lock.write_all(body.as_bytes()) {
                eprintln!("write stdout: {}", e);
                return ExitCode::from(1);
            }
        }
    }
    ExitCode::SUCCESS
}

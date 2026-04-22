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
  --from PATH     resume from a previously-emitted JSON state file; extends
                  the saved computation up to MAX_STEM.  MAX_STEM must be
                  ≥ the saved max_stem.  Output format / destination still
                  controlled by the flags above.
  -o, --output PATH
                  write output to PATH instead of stdout.  Turns off ANSI
                  unless --color is also given.
  -h, --help      show this help.

examples:
  table 12                                 # colored report to terminal
  table 12 --plain -o report.txt           # save .txt report
  table 12 --ascii -o report.txt           # strict ASCII .txt report
  table 12 --json -o report.json           # JSON for Python visualizer
  table 26 --from state24.json --json -o state26.json
                                           # resume 24 → 26
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
    resume_from: Option<PathBuf>,
}

fn parse_args() -> Result<Args, String> {
    let mut max_stem: Option<usize> = None;
    let mut format = Format::Human;
    let mut style: Option<RenderStyle> = None;
    let mut output: Option<PathBuf> = None;
    let mut resume_from: Option<PathBuf> = None;

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
            "--from" => {
                i += 1;
                if i >= argv.len() {
                    return Err(format!("{} needs an argument", a));
                }
                resume_from = Some(PathBuf::from(&argv[i]));
            }
            s if s.starts_with("--from=") => {
                resume_from = Some(PathBuf::from(&s["--from=".len()..]));
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
        resume_from,
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

    let table = match &args.resume_from {
        Some(path) => {
            let src = match std::fs::read_to_string(path) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("read {}: {}", path.display(), e);
                    return ExitCode::from(1);
                }
            };
            let mut t = match CurtisTable::from_json(&src) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("parse {}: {}", path.display(), e);
                    return ExitCode::from(1);
                }
            };
            if args.max_stem < t.max_stem {
                eprintln!(
                    "MAX_STEM ({}) < saved max_stem ({}); cannot shrink.",
                    args.max_stem, t.max_stem
                );
                return ExitCode::from(2);
            }
            if args.max_stem == t.max_stem {
                eprintln!(
                    "Loaded state from {} (stem {}); nothing to extend.",
                    path.display(), t.max_stem
                );
            } else {
                eprintln!(
                    "Loaded state from {} (stem {}); extending to stem {}...",
                    path.display(), t.max_stem, args.max_stem
                );
                t.extend_to(args.max_stem);
            }
            t
        }
        None => {
            eprintln!("Computing Curtis table through stem {}...", args.max_stem);
            CurtisTable::compute(args.max_stem)
        }
    };

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

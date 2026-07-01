mod cli;

use clap::Parser;

fn main() -> std::process::ExitCode {
    cli::Cli::parse().run()
}

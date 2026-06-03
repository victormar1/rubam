//! Crate error type. Public so HARMOS can match on it.

use thiserror::Error;

/// Result alias for the rubam public API.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors produced by the rubam public API.
#[derive(Debug, Error)]
pub enum Error {
    /// Underlying I/O error (open, read, seek).
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// The reference name was not in the BAM header.
    #[error("reference {0:?} not in BAM header")]
    UnknownReference(String),

    /// A CIGAR operation could not be decoded.
    #[error("CIGAR decode error: {0}")]
    Cigar(String),

    /// An auxiliary tag could not be decoded.
    #[error("aux decode error: {0}")]
    Aux(String),
}

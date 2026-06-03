//! Header — TID-indexed access to BAM reference sequences.

use noodles::sam;

/// Owned wrapper around `noodles::sam::Header` with HARMOS-shaped helpers.
#[derive(Clone, Debug)]
pub struct Header {
    inner: sam::Header,
}

impl Header {
    /// Build a `Header` from a parsed `noodles::sam::Header`.
    pub fn from_noodles(inner: sam::Header) -> Self {
        Self { inner }
    }

    /// The underlying `noodles::sam::Header`. Useful when calling deeper APIs.
    pub fn as_noodles(&self) -> &sam::Header {
        &self.inner
    }

    /// Number of reference sequences. Returns 0 if the header has no `@SQ` lines.
    pub fn target_count(&self) -> u32 {
        self.inner.reference_sequences().len() as u32
    }

    /// Reference name for `tid`. Returns `None` if `tid` is out of bounds.
    pub fn tid2name(&self, tid: u32) -> Option<&[u8]> {
        self.inner
            .reference_sequences()
            .get_index(tid as usize)
            .map(|(name, _)| name.as_ref())
    }

    /// Reference length (in bp) for `tid`. Returns `None` if `tid` is out of bounds.
    pub fn target_len(&self, tid: u32) -> Option<u64> {
        self.inner
            .reference_sequences()
            .get_index(tid as usize)
            .map(|(_, ref_seq)| ref_seq.length().get() as u64)
    }

    /// Iterator over reference names in the order they appear in `@SQ` lines.
    pub fn target_names(&self) -> impl Iterator<Item = &[u8]> + '_ {
        self.inner
            .reference_sequences()
            .iter()
            .map(|(name, _)| name.as_ref())
    }
}

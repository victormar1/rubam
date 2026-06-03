// tests/api_aux.rs
use rubam::api::Aux;

#[test]
fn variant_constructors_compile() {
    let _ = Aux::Char(b'A');
    let _ = Aux::I8(-1);
    let _ = Aux::U8(255);
    let _ = Aux::I16(-30000);
    let _ = Aux::U16(60000);
    let _ = Aux::I32(-1_000_000);
    let _ = Aux::U32(3_000_000_000);
    let _ = Aux::Float(0.5);
    let _ = Aux::Double(1.5);
    let _ = Aux::String("hello");
    let _ = Aux::HexByteArray("DEADBEEF");
    let _ = Aux::ArrayI8(&[-1, 0, 1]);
    let _ = Aux::ArrayU8(&[0, 1, 255]);
    let _ = Aux::ArrayI16(&[-1, 0, 1]);
    let _ = Aux::ArrayU16(&[0, 1, 65000]);
    let _ = Aux::ArrayI32(&[-1, 0, 1]);
    let _ = Aux::ArrayU32(&[0, 1, 4_000_000_000]);
    let _ = Aux::ArrayFloat(&[1.0, 2.0]);
}

#[test]
fn aux_string_pattern_match() {
    // The HARMOS use case: pattern-match on Aux::String
    let a: Aux = Aux::String("chr21,1000,+,50M,60,0;");
    if let Aux::String(s) = a {
        assert!(s.starts_with("chr21,"));
    } else {
        panic!("expected Aux::String");
    }
}

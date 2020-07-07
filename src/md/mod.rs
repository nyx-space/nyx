pub mod ui;

/// A Mission Design handler
pub trait MdHdlr<StateType: Copy> {
    fn handle(&mut self, state: &StateType);
}

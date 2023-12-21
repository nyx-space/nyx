use egui;
mod app;
mod controls;

pub use app::NyxGui;

pub trait View {
    fn ui(&mut self, ui: &mut egui::Ui);
}

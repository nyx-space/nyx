/// This module defines the top-level UI app that will be called at main.rs
use eframe::{App, CreationContext};

use crate::gui::controls::ScenarioPicker;

use super::View;

pub struct NyxGui;

impl Default for NyxGui {
    fn default() -> Self {
        Self {
            // specify defaults here
        }
    }
}

impl NyxGui {
    pub fn new(cc: &CreationContext<'_>) -> Self {
        // This is also where you can customize the look and feel of egui using
        // `cc.egui_ctx.set_visuals` and `cc.egui_ctx.set_fonts`.

        Default::default()
    }
}

impl App for NyxGui {
    /// By default the app is Reactive: egui is only updated if are input events
    /// (like mouse movements) or there are some animations in the GUI.
    fn update(&mut self, ctx: &egui::Context, frame: &mut eframe::Frame) {
        egui::CentralPanel::default().show(ctx, |ui| {
            ScenarioPicker::default().ui(ui);
        });
    }
}

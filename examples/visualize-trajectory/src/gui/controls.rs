use crate::scenario::{lunar_transfer, Scenarios};
pub struct ScenarioPicker {
    scenario: Scenarios,
}

impl Default for ScenarioPicker {
    fn default() -> Self {
        Self {
            scenario: Scenarios::LunarTransfer,
        }
    }
}

impl super::View for ScenarioPicker {
    fn ui(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            // run the scenario when button is clicked
            if ui.button("go").clicked() {
                match self.scenario {
                    Scenarios::LunarTransfer => lunar_transfer(),
                    _ => {}
                }
            };
        });
    }
}

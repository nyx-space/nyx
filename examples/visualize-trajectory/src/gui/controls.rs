use crate::scenario::Scenarios;
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
            // change the scenario enum from a dropdown
            egui::ComboBox::from_label("Scenario")
                    .selected_text(self.scenario.to_string())
                    .show_ui(ui, |ui| {
                        for scen in &[
                            Scenarios::LunarTransfer,
                            Scenarios::OrbitDesign,
                        ] {
                            ui.selectable_value(&mut self.scenario, *scen, scen.to_string());
                        }
                    });

            // run the scenario when button is clicked
            if ui.button("run").clicked() {
                self.scenario.run()
            };
        });
    }
}

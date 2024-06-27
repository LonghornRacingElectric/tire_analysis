from tire_model.file_processing._process_tir import _Processor

new_tire = _Processor(name = "test", file_path = "./tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_8psi_PAC02_UM2.tir")

params = new_tire._import_data(path = "./tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_8psi_PAC02_UM2.tir")

print(params)
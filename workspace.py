from tire_toolkit import workflows

tires = workflows.import_tires(
    tire_names = 
        [
        "Hoosier_16_7.5x7_R25B_8psi",
        "Hoosier_16_7.5x7_R25B_10psi",
        "Hoosier_16_7.5x7_R25B_12psi",
        "Hoosier_16_7.5x7_R25B_14psi"
        ],
    tir_paths = 
        [
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_8psi_PAC02_UM2.tir",
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_10psi_PAC02_UM2.tir",
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_12psi_PAC02_UM2.tir",
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_14psi_PAC02_UM2.tir"
        ]
)

workflows.plot_tires(
    tire_object = tires, 
    names = [
        "Hoosier_16_7.5x7_R25B_8psi",
        "Hoosier_16_7.5x7_R25B_10psi",
        "Hoosier_16_7.5x7_R25B_12psi",
        "Hoosier_16_7.5x7_R25B_14psi"], 
    FZ_min = 100, 
    FZ_max = 1091, 
    file_name = "First_Principles_Metrics.pdf")
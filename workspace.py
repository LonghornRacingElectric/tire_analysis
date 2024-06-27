from tire_toolkit import workflows

tires = workflows.import_tires(
    tire_names = 
        [
        "Hoosier_16_7.5x7_LC0"
        ],
    tir_paths = 
        [
        "./tire_toolkit/assets/tir_files/Round_8_Hoosier_R25B_16x7p5_10_on_7in_8psi_PAC02_UM2.tir"
        ]
)

fig_1 = tires.single_plot(name = "Hoosier_16_7.5x7_LC0", FZ_min = 100, FZ_max = 1091)

workflows.save_pdf(figures = [fig_1], file_name = "First_Principles_Metrics.pdf", save_location = "./outputs")
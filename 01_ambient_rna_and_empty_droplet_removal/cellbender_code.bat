::CellBender background removal for scRNA-seq datasets

:: naive - batch 1
cellbender remove-background --input "\path\to\Naive_B1\raw_feature_bc_matrix.h5" --output "\path\to\Naive_B1\Naive_B1" 

:: cci_chronic - batch 3
cellbender remove-background --input "\path\to\CCI_Chronic_B3\raw_feature_bc_matrix.h5" --output "\path\to\CCI_Chronic_B3\CCI_Chronic_B3" 

:: cci_mid - batch 3
cellbender remove-background --input "\path\to\CCI_Mid_B3\raw_feature_bc_matrix.h5" --output "\path\to\CCI_Mid_B3\CCI_Mid_B3"

:: cci_chronic - batch 6
cellbender remove-background --input "\path\to\CCI_Chronic_B6\raw_feature_bc_matrix.h5" --output "\path\to\CCI_Chronic_B6\CCI_Chronic_B6" 

:: sni_acute - batch 6
cellbender remove-background --input "\path\to\SNI_Acute_B6\raw_feature_bc_matrix.h5" --output "\path\to\SNI_Acute_B6\SNI_Acute_B6" 



# Glycerol Estimation and Numerical Analysis Using MATLAB

## Project Overview

This project focuses on glycerol concentration estimation using MATLAB-based numerical analysis and statistical modeling techniques. Multiple experimental datasets are processed and analyzed to investigate the relationship between impedance-derived features and glycerol concentration. Model outputs are compared against reference measurements from the original dataset to evaluate calibration accuracy and robustness.

The study applies preprocessing, signal alignment, feature extraction, and multivariate statistical methods to experimental QCM flow data.

---

## Objectives

* Analyze QCM impedance and flow-based experimental datasets
* Process and align multiple dataset sources (`data2`, `data3`, `data4`, `IWIS2025flow`)
* Estimate glycerol concentration using regression-based models
* Compare estimated values with reference measurements
* Evaluate calibration performance using statistical metrics
* Visualize relationships between features and concentration

---

## Datasets

The project utilizes the following datasets:

| Dataset        | Purpose                                                                |
| -------------- | ---------------------------------------------------------------------- |
| `data2`        | Processed experimental dataset                                         |
| `data3`        | Processed experimental dataset with additional MATLAB analysis scripts |
| `data4`        | Processed experimental dataset                                         |
| `IWIS2025flow` | Reference dataset containing flow and impedance measurements           |

These datasets consist of impedance-derived features, flow control signals, and auxiliary measurement parameters related to glycerol-in-water experiments.

---

## Methodology

### Data Preprocessing

The preprocessing pipeline includes:

* Time alignment between flow and impedance measurements
* Removal of corrupted or unstable measurement regions
* Feature extraction from impedance signals
* Standardization and numerical normalization
* Organization of steady-state measurement windows

---

### Statistical and Numerical Analysis

* Segmentation of experimental data into steady-state regions
* Computation of mean and standard deviation per segment
* Construction of calibration-ready feature matrices

---

### Regression and Multivariate Modeling

Regression techniques are used to estimate glycerol concentration from impedance-derived features:

* Linear Regression (univariate calibration)
* Principal Component Analysis (PCA) for dimensionality reduction
* Partial Least Squares (PLS) regression for multivariate calibration
* Model comparison between univariate and multivariate approaches

---

### Model Evaluation

Model performance is evaluated using:

* Calibration sensitivity (Hz/%)
* Coefficient of determination (R²)
* Error metrics (MAE, MSE, RMSE)
* Noise floor estimation
* Limit of detection (LOD) analysis

---

## Technologies Used

| Technology                              | Purpose                                                                       |
| --------------------------------------- | ----------------------------------------------------------------------------- |
| MATLAB                                  | Numerical computation, signal processing, statistical modeling, visualization |
| Statistics and Machine Learning Toolbox | Regression and multivariate analysis                                          |
| Signal Processing Toolbox               | Feature extraction and signal analysis                                        |

---

## Project Structure

```text
numerical_project/
│
├── IWIS2025flow/
│   └── data.txt
│
├── data2/
│   ├── B1_data.txt
│   ├── B2_data.txt
│   ├── G_data.txt
│   ├── R_data.txt
│   ├── X1_data.txt
│   ├── X2_data.txt
│   ├── Y_abs2_data.txt
│   ├── Y_abs_data.txt
│   ├── Z_abs2_data.txt
│   ├── Z_abs_data.txt
│   ├── angle_Z_data.txt
│   └── pp.mat
│
├── data3/
│   ├── B1_data.txt
│   ├── B2_data.txt
│   ├── G_data.txt
│   ├── LAG.m
│   ├── PCA.m
│   ├── R_data.txt
│   ├── X1_data.txt
│   ├── X2_data.txt
│   ├── Y_abs2_data.txt
│   ├── Y_abs_data.txt
│   ├── Z_abs2_data.txt
│   ├── Z_abs_data.txt
│   ├── angle_Z_data.txt
│   ├── mdpi_3rd_poly.m
│   ├── mdpi_poly_data.m
│   ├── pls.m
│   ├── pp.mat
│   └── pp_withflow.mat
│
└── data4/
    ├── B1_data.txt
    ├── B2_data.txt
    ├── G_data.txt
    ├── R_data.txt
    ├── X1_data.txt
    ├── X2_data.txt
    ├── Y_abs2_data.txt
    ├── Y_abs_data.txt
    ├── Z_abs2_data.txt
    ├── Z_abs_data.txt
    ├── angle_Z_data.txt
    └── pp.mat
```

---

## Requirements

The project is implemented entirely in MATLAB.

Required software:

* MATLAB R2021a or newer
* Statistics and Machine Learning Toolbox
* Signal Processing Toolbox

---

## Usage

Open MATLAB and navigate to the project directory:

```matlab
numerical_project/IWIS2025flow
```

Run analysis scripts:

```matlab
run('mdpi_poly_data.m')
run('PCA.m')
run('pls.m')
```

---

## Evaluation Metrics

* Mean Absolute Error (MAE)
* Mean Squared Error (MSE)
* Root Mean Squared Error (RMSE)
* R² Score

These metrics quantify agreement between estimated and reference glycerol values.

---

## Visualization Outputs

* Predicted vs reference comparison plots
* Correlation heatmaps
* Residual analysis
* Regression calibration curves
* Feature distribution plots

---

## Results

The analysis demonstrates a consistent relationship between impedance-derived features and glycerol concentration. Multivariate approaches improve robustness against feature collinearity compared to univariate calibration models. Overall performance is evaluated through calibration accuracy, error reduction, and stability across experimental conditions.

---

## License

This project is intended for academic, research, and educational use.

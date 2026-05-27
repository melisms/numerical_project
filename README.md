# Glycerol Prediction and Numerical Analysis Using Machine Learning

## Project Overview

This project focuses on predicting glycerol values through machine learning and numerical analysis techniques. Multiple processed datasets are used to train regression models, and the generated predictions are compared against the original dataset to evaluate model performance and prediction accuracy.

The study investigates the relationship between transformed datasets and the original glycerol measurements while applying statistical analysis, preprocessing, visualization, and supervised learning methods.

---

## Objectives

* Analyze glycerol-related numerical datasets
* Train machine learning regression models using processed data
* Predict glycerol values from experimental features
* Compare predicted results with original measurements
* Evaluate model accuracy using regression metrics
* Visualize relationships and prediction performance

---

## Datasets

The project utilizes multiple datasets:

| Dataset         | Purpose                          |
| --------------- | -------------------------------- |
| `data2`         | Training dataset                 |
| `data3`         | Training dataset                 |
| `data4`         | Training dataset                 |
| `IWIS2025flow`  | Reference dataset for comparison |

These datasets contain numerical measurements associated with glycerol values and experimental parameters.

---

## Methodology

### Data Preprocessing

The preprocessing pipeline includes:

* Missing value handling
* Data cleaning
* Feature extraction
* Normalization and scaling
* Numerical formatting

### Machine Learning Workflow

The machine learning process consists of:

1. Loading processed datasets
2. Splitting training and testing data
3. Training regression models
4. Generating glycerol predictions
5. Comparing predictions with original values
6. Evaluating model performance

### Regression Analysis

Regression techniques are applied to estimate glycerol concentrations and numerical trends from processed experimental data.

Potential models include:

* Linear Regression
* Random Forest Regression
* Gradient Boosting Regression
* Support Vector Regression
* Ensemble-based approaches

---

## Technologies Used

| Technology | Purpose                                                                                          |
| ---------- | ------------------------------------------------------------------------------------------------ |
| MATLAB     | Numerical analysis, machine learning, data preprocessing, regression modeling, and visualization |

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

The project is developed and executed entirely in MATLAB.

Required software:

* MATLAB R2021a or newer
* Statistics and Machine Learning Toolbox
* Signal Processing Toolbox

---

## Usage

1. Open MATLAB.

2. Clone or download the repository:

```bash
git clone https://github.com/melisms/numerical_project.git
```

3. Open the project folder in MATLAB:

```text
numerical_project/
```

4. Navigate to the `IWIS2025flow` directory.

5. Run the required MATLAB scripts for preprocessing, regression analysis, machine learning training, and prediction.

Example:

```matlab
run('mdpi_poly_data.m')
run('PCA.m')
run('pls.m')
```

6. Compare predicted glycerol values with the original dataset outputs.

---

## Evaluation Metrics

The project evaluates prediction performance using:

* Mean Absolute Error (MAE)
* Mean Squared Error (MSE)
* Root Mean Squared Error (RMSE)
* R² Score

These metrics are used to measure how closely predicted glycerol values match the original experimental measurements.

---

## Visualization Outputs

The analysis may generate:

* Predicted vs Actual plots
* Correlation heatmaps
* Residual analysis plots
* Regression comparison graphs
* Distribution visualizations

---

## Results

The trained models are evaluated by comparing predicted glycerol values against the original dataset. Lower error metrics and higher R² values indicate improved predictive performance and stronger agreement with the reference measurements.

---

## License

This project is intended for academic, research, and educational purposes.

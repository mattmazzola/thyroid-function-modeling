## Thyroid Function Modeling (Streamlit)

This Streamlit app calculates SPINA-GD, SPINA-GT, and Jostel's TSH Index (TSHI) from user-entered lab values with selectable units. Constants and molecular weights are adjustable.

### Features
- Inputs with unit pickers for Free T3, Free T4, Reverse T3, TSH
- Adjustable molecular weights (T3, T4) and model constants
- Outputs: SPINA-GD (nmol/s), SPINA-GT (pmol/s), TSHI
- Help expanders for each metric and interpretation hints

### Formulas
- SPINA-GT (using FT4):
  GT = βT × (DT + [TSH])/(αT × [TSH]) × (1 + K41[TBG] + K42[TBPA]) × [FT4]
  (pmol/s when converted)

- SPINA-GD (using FT3, FT4):
  GD = β31 × (KM1 + [FT4])/(1 + K30[TBG]) × [FT3]/(α31 × [FT4])
  (nmol/s when converted)

- Jostel's TSH Index (TSHI):
  TSHI = ln(TSH[mIU/L]) + 0.1345 × FT4[pmol/L]

References: Dietrich et al., Front Endocrinol (Lausanne). 2016; SPINA Thyr documentation.

### Quick start

```bash
pip install -r requirements.txt
streamlit run app.py
```

### Notes
Results are for educational purposes and should not replace clinical judgment.


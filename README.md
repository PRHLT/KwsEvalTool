# KwsEvalTool

Compute *(interpolated-)precision and recall curve* (**R-P**), *average precision* (**AP**) and/or *mean average precision* (**mAP**) for evaluating keyword spotting performance.

## Compilation

```bash
gcc -Wall -O4 -o kws-assessment kws-assessment.c
```

## Usage

To compute **AP** and **mAP**:
```bash
kws-assessment -t -s -a -m data_test.dat -w egs/keywords_test.lst \
               -l 16376 egs/data_test.dat
```

To generate the corresponding **R-P** curve:
```bash
kws-assessment -t -s data_test.dat -w egs/keywords_test.lst \
               -l 16376 egs/data_test.dat > egs/r-p_data.dat
```
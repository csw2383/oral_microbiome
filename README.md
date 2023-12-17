# Comparative Analysis of Human Oral Microbiome Communities in Health and Disease Project

There are 20 samples we gather from different invidual's cavity with different health condition which ends with the suffix

```
.motus.ratio.genus
```

The code and figures in the repo are mainly for the result and analysis part of the final project. Running the following code in this directory will allow you to generate the pie chart for your data. However, you do have to specify your path in the code of where your data set (that ends with the motus.ratio.genus) is. After that, all figures will be automatically generated showing the relative frequency of all bacteria genera in the sample. There is also a CSV name sample_name which defines which sample we are looking at and the condition of the sample. 



```
python main.py
```

Make sure to modify the directory to yours:
```
directory = '/Users/weiwei/cu/cs4775/motus_output'  # Replace this with your directory path
```

The Code also includes the PCA where running the command

```
python PCA_plot_genus.py
```

will generate several PCA figures which is something we work on in the workdir at BIOHPC



### 1) Starting from a list of samples:

```
GM12878"
B cell (CD19-positive)
```

### 2) Select TSS from samples

```
from GUD import GUDglobals
from GUD.scripts.sample2gene import get_differentially_expressed_tss

# Establish a GUD session
session = GUDglobals.establish_GUD_session()

# Get differentially expressed TSSs
diff_exp_tss = get_differentially_expressed_tss(
        session,
        samples=["GM12878", "B cell (CD19-positive)"]
    )
```    


3) Select a region for a gene:

```
from GUD import GUDglobals
from GUD.scripts.gene2region import get_gene_region

# Establish a GUD session
session = GUDglobals.establish_GUD_session()

# Get that gene's region
region = get_gene_region(
    session,
    gene,
    samples=["GM12878", "B cell (CD19-positive)"]
)
```
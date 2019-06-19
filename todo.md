## Genomic Features 

### clinvar

```
## queries 
select_by_location(cls, session, chrom, start, end, as_genomic_feature=False)
select_by_clinvarID(cls, session, clinvarID, as_genomic_feature=False)
is_unique(cls, session, clinvarID)
```

- [ ] lift select_by_location 
- [ ] add/lift select_exact_location
- [ ] lift is_unique
- [ ] add/lift select_by_uid
- [ ] remove option for as_genomic_feature, make automatic

### conservation

```
## queries 
select_by_location(cls, session, chrom,start, end, as_genomic_feature=False)
is_unique(cls, session, regionID,sourceID) 
select_unique(cls, session, regionID, sourceID)
```

- [ ] lift select_by_location 
- [ ] add/lift select_exact_location
- [ ] lift is_unique
- [ ] add/lift select_by_uid
- [ ] remove option for as_genomic_feature, make automatic

### copy_number_variants

```
## queries 
select_by_location(cls, session, chrom, start, end, as_genomic_feature=False)
select_by_exact_location(cls, session, chrom, start, end,as_genomic_feature=False)
select_by_uid(cls, session, uid, as_genomic_feature=False)
is_unique(cls, session, name) 
```

- [ ] lift select_by_location 
- [ ] lift select_exact_location
- [ ] lift is_unique
- [ ] lift select_by_uid
- [ ] remove option for as_genomic_feature, make automatic

### dna_accessibility

```
## queries 
is_unique(cls, session, regionID, sampleID, experimentID, sourceID)
select_unique(cls, session, regionID, sampleID, experimentID, sourceID)
```

- [ ] make genomic feature

### enhancers

```
## queries 
is_unique(cls, session, regionID, sampleID,experimentID, sourceID)
select_unique(cls, session, regionID,sampleID, experimentID, sourceID)
select_by_sample(cls, session,sample, as_genomic_feature=False)
```

- [ ] lift is_unique
- [ ] lift select_unique
- [ ] remove as_genomic_feature option 

### genes

```
## queries 
is_unique(cls, session, regionID, name,sourceID)
select_unique(cls, session, regionID,name, sourceID)
select_by_location(cls, session, chrom,start, end, as_genomic_feature=False)
select_by_name(cls, session, name,as_genomic_feature=False)
select_by_names(cls, session, names=[],as_genomic_feature=False)
select_by_uid(cls, session, uid,as_genomic_feature=False)
```

- [ ] remove as_genomic_feature option 
- [ ] remove select by name 
- [ ] lift select_by_uid
- [ ] lift select_by_location 
- [ ] lift select_unique
- [ ] lift is_unique 

### histone_modifications

```
## queries 
select_by_location(cls, session, chrom, start, end)
is_unique(cls, session, regionID, sourceID)
```

- [ ] make into genomic feature
- [ ] lift select_by_location
- [ ] lift is_unique 

### rmsk

```
## queries 
```

- [ ] 

### short_tandem_repeats

```
## queries 
```

- [ ] 

### tads

```
## queries 
```

- [ ] 

### tf_binding

```
## queries 
```

- [ ] 

### transcription_start_sites

```
## queries 
```

- [ ] 

## Non-Genomic Features 

### chroms

```
## queries 
select_by_chrom(cls, session, chrom)        ##get rid of 
select_by_chroms(cls, session, chroms=[])
chrom_size(cls, session, chrom)             ## get rid of 
chrom_sizes(cls, session, chroms=[])

```

- [ ] get rid of select by chrom and chrom_size

### experiments

```
## queries 
```

- [ ] 

### expression

```
## queries 
```

- [ ] 

### regions

```
## queries 
```

- [ ] 

### samples

```
## queries 
```

- [ ] 

### sources

```
## queries 
```

- [ ] 


# GUD API

## Routes

| Method | Endpoint                                                     | Description                                               |
| :----- | ------------------------------------------------------------ | --------------------------------------------------------- |
| GET    | /resource?chrom={str}&start={int}&end={int}&location={str}   | Return a list resources* by location                      |
| GET    | /resource?chrom={str}&start={int}&end={int}&location={str}&sources={str} | Return a list resources* by location and source           |
| GET    | /resource?chrom={str}&start={int}&end={int}&location={str}&samples={str} | Return a list of resources^ by location and sample        |
| GET    | /resource?chrom={str}&start={int}&end={int}&location={str}&experiments={str} | Return a list of resources** by location and experiment   |
| GET    | /genes?names={str}                                           | Return a list of genes by gene symbol                     |
| GET    | /genesymbols                                                 | Return all gene symbols                                   |
| GET    | /clinvar?clinvar_id={str}                                    | Return clinvar variant by clinvarID                       |
| GET    | /short_tandem_repeats/pathogenic                             | Return all pathogenic short tandem repeats                |
| GET    | /short_tandem_repeats?motif={str}&rotation={bool}            | Return all short tandem repeats matching a specific motif |
| GET    | /tads?chrom={str}&start={int}&end={int}&location={str}&restriction_enzymes={str} | Return all tads by location and restriction enzyme        |
| GET    | /tf_binding?chrom={str}&start={int}&end={int}&location={str}&tfs={str} | Return all tf binding sites by location and tf            |
| GET    | /tss?genes={str}                                             | Return all tss by genes                                   |
| GET    | /genic_tss                                                   | Return all genic tss                                      |
| GET    | /sources                                                     | return all sources                                        |
| GET    | /samples                                                     | Return all samples                                        |
| GET    | /experiments                                                 | Return all experiments                                    |
| GET    | /chroms                                                      | Return all chromosomes                                    |


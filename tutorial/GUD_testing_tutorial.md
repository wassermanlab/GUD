# GUD - Genomic Unification Database Testing Tutorial

### Background

The Genomic Unification Database (GUD) is a project developed by the Wasserman lab to house a unified version of several different data-sets and data-types describing human genomic information. The database is designed to reduce data redundancies found in the original data-sets and to leverage the full capacities of mySQL. The back end code is in Python and Flask and uses a mySQL database housed on the CMMT servers. For ease of access we have developed a simple API for querying the database which should be used to retrieve data.   

![](/home/tamar/Desktop/GUD/tutorial/pics/Fig.2.png)

### Pre-requisites 

- `CMMT VPN`

- browser 
- programming language of choice for api access

### Available Resources 

GUD follows this basic schema, from the database you can fetch all the Genomic Feature 1 and 2 (GF1/GF2) features, chroms, samples, experiments, and sources. For the purpose of this tutorial we have inserted a subset of the data, mainly these tables: `chroms, clinvar, conservation, copy_number_variants, cpg_islands, experiments, genes, histone_modifications, regions, rmsk, samples, short_tandem_repeats, sources, tf_binding`.

![](/home/tamar/Desktop/GUD/tutorial/pics/GUD_schema-orm.png)

### Access Model 

GUD is designed for programmatic access through an API. A front end API website is available at the URL `http://gud.cmmt.ubc.ca:8080/`, for now this is only accessible if you are on the CMMT network, either through a VPN or via a Ethernet cable. This website has some basic documentation on access as well as a live API where users can test small queries to understand the query syntax. Please use the website for learning how to use the API. 

![](/home/tamar/Desktop/GUD/tutorial/pics/Screenshot from 2019-11-25 16-58-41.png)

### API Access Scripts

For faster programmatic access it is best to fetch blocks of data through scripts. We provide simple scripts for querying the database programmatically in `R` and `Python`. For other languages please consult google! Currently the API caps request from the same IP address at 5 per second so please be aware when designing your own scripts. 

#### R

These are scripts for loading results into a data frame. Keep in mind this might take a while if you are trying to fetch millions of rows.

```R
# if not installed uncomment and install these packages
# install.packages("httr")
# install.packages("jsonlite")

# requirements
require("httr")
require("jsonlite")

# fetching a single page
base <- "http://gud.cmmt.ubc.ca:8080/api/v1/"
db <- "hg38/"
resource <- "genes" # replace this with your desired query resource and filters
url <- paste(base, db, resource, sep="")
page <- GET(url)
page_text <- content(page, "text")
page_json <- fromJSON(page_text, flatten = TRUE)

# fetching all pages from a query 
full_set <- data.frame() 
base <- "http://127.0.0.1:5000/api/v1/"
db <- "hg38/"
resource <- "short_tandem_repeats?pathogenicity=true"   # replace this with your desired query resource and filters
# add whatever filters you want to resource   
next_url <- paste(base, db, resource, sep="")

while(!is.null(next_url)){
  page <-GET(next_url)
  page_text <- content(page, "text", encoding="UTF-8")
  page_json <- fromJSON(page_text, flatten = TRUE)
  next_url <- page_json$`next`
  results <- page_json$results
  full_set <- rbind(full_set, results)
}

```

#### Python

Requirements: 

- Python 3
- requests module, install with`pip install requests`

The following script takes an initial script and gets all the results in an array, this can later be parsed into whatever data structure the user wants. Keep in mind this might take a while if you are trying to fetch millions of rows.

```python
import requests
import json
import time

def get_all_results(request_url):
    results = []
    next = True
    url = request_url
    # in this while loop wait 1 sec
    counter = 0 
    while(next):
        if (counter == 4):
            counter = 0 
            time.sleep(1)
        else: 
            counter = counter + 1
        req = requests.get(url)
        if req.status_code == 200:
            res = json.loads(req.content.decode('utf-8'))
            results = results + res['results']
            if ('next' in res):
                url = res['next']
            else: 
                next = False
        else:
            return False
    return results
```

### Testing Guide

#### Steps: 

1. Check that you are on the CMMT network, either via wired connection or VPN

2. Navigate to http://gud.cmmt.ubc.ca:8080/

3. Test navigation bar make sure you can navigate to all places)

   1. http://gud.cmmt.ubc.ca:8080/contact
   2. http://gud.cmmt.ubc.ca:8080/docs
   3. http://gud.cmmt.ubc.ca:8080/live_api

4. In http://gud.cmmt.ubc.ca:8080/docs test that the filtering works

5. In http://gud.cmmt.ubc.ca:8080/docs make sure everything is understandable 

6. Test http://gud.cmmt.ubc.ca:8080/live_api (expect results only from the following tables until further notice: genes, chroms, samples, short_tandem_repeats, copy_number_variants, clinvar)

    **NOTE:** Be aware that GUD expects 1 based inputs and returns 0 based coordinates. 

#### Bugs:

- We are expecting bugs to arise when you test the live_api, help us by reporting all bugs here: https://docs.google.com/spreadsheets/d/1LxgiacPJcwV7YU5gEHJFkRBdj7Ac_XxlRR-u1GycnZk/edit?usp=sharing. 
- When reporting bugs please write the exact steps to reproduce the bug
- Include screenshots if the bug is something visual 
- Optional but helpful - if you can open the devtools (ctrl-shift-i on chrome) and check the console tab for any errors and record them

#### Some examples of bugs:

- Sending a request and getting “loading...” for more than a 5 minutes

- Making an invalid query and getting back a useless message i.e. `"error": "The server encountered an internal error and was unable to complete your request. Either the server is overloaded or there is an error in the application."`

- Getting a response that is wrong i.e. if you search for a gene and the coordinates are clearly wrong

- Any visual problems 

#### Comments: 

Any general comments/suggestions are highly appreciated, please message me (Tamar) on slack or by email tavshalom@cmmt.ubc.ca. 
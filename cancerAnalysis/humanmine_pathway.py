#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The following two lines will be needed in every python script:
from intermine.webservice import Service
import pandas as pd

service = Service("https://www.humanmine.org/humanmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "primaryIdentifier", "symbol", "proteins.primaryAccession",
    "proteins.primaryIdentifier", "proteins.pathways.identifier",
    "proteins.pathways.name", "proteins.pathways.dataSets.name"
)

# Uncomment and edit the line below (the default) to select a custom sort order:
query.add_sort_order("Gene.primaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("proteins.pathways.name", "IS NOT NULL", code="A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A")

allData = {}
allData["primaryIdentifier"] = []
allData["symbol"] = []
allData["proteins.primaryAccession"] = []
allData["proteins.primaryIdentifier"] = []
allData["proteins.pathways.identifier"] = []
allData["proteins.pathways.name"] = []
allData["proteins.pathways.dataSets.name"] = []

for row in query.rows():
    allData["primaryIdentifier"].append(row["primaryIdentifier"])
    allData["symbol"].append(row["symbol"])
    allData["proteins.primaryAccession"].append(row["proteins.primaryAccession"])
    allData["proteins.primaryIdentifier"].append(row["proteins.primaryIdentifier"])
    allData["proteins.pathways.identifier"].append(row["proteins.pathways.identifier"])
    allData["proteins.pathways.name"].append(row["proteins.pathways.name"])
    allData["proteins.pathways.dataSets.name"].append(row["proteins.pathways.dataSets.name"])

allData_df = pd.DataFrame(allData)
print(allData_df)
allData_df.to_csv("humanmine_pathway_AllData.csv", index=False)

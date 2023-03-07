#############################################################################
#### The following script downloads files from several public databases.
#### WARNING: some of the databases are protected by the electronic user agreement. By downloading the files, the user agrees to be bound to the terms and conditions of the corresponding electronic agreement. Checking the terms and conditions of each agreement is the responsibility of the end user not the owners of the following script.
### In case of the databases update, the old urls might stop working. You will need to obtain the most recent version of the database.


# REGULONDB 11.1

url = "https://regulondb.ccg.unam.mx/menu/download/datasets/files/NetWorkTFGene.txt" # in case url stops working update the link by locating TF - gene interactions dataset at Downloads/Experimental Datasets tab at https://regulondb.ccg.unam.mx/index.jsp

# The file is downloaded from public database RegulonDB 11.1 under an Academic/Noncommercial Use License. By downloading the file below, the user agrees to be bound to the terms and conditions of the electronic End User License Agreement, you can find the corresponding agreement here: https://regulondb.ccg.unam.mx/menu/download/full_version/terms_and_conditions.jsp


destfile = 'data_source/network_tf_gene.txt'  
download.file(url, destfile)


# SUBTIWIKI version from 2023-02-27

url = "http://subtiwiki.uni-goettingen.de/v4/regulation/exporter" # in case url stops working update the link by locating Regulations dataset at Exports tab at http://subtiwiki.uni-goettingen.de/

destfile = 'data_source/Subtiwiki.csv'  
download.file(url, destfile)



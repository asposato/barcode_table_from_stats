from collections import Counter

stats_file = input("Type the name of your stats file including the extension (i.e. test.stats) and press Enter/Return: ")
stats_summary_file_almost = stats_file.split(".")
stats_file_summary = str(stats_summary_file_almost[0]) + '_IDs_umi_barcodes.tsv'
output_all_barcodes = str(stats_summary_file_almost[0]) + '_all_barcodes.tsv'
output_top_barcodes = str(stats_summary_file_almost[0]) + '_top_barcode.tsv'

print('Writing intermediate summary table of stats file...')
output = open(stats_file_summary, 'w')
output.write('cell.id' + '\t' + 'umi' + '\t' + 'PASS' + '\t'+'site1' + '\t' + 'site2'+ '\t' + 'site3' + '\t' + 'site4' + '\t' + 'site5' + '\t' + 'site6' + '\t' + 'site7' '\t' + 'site8' + '\t' + 'site9' + '\t' + 'site10' + '\t'+'\n')

stats = open(stats_file, "r")
stats.readline()
for line in stats.readlines():
    read = line.split("_")
    cell_id = read[1][:16]
    umi = read[1][16:26]
    barcode = read[6]
    barcode = barcode.split("\t")
    PASS = barcode[1]
    site1 = barcode[22]
    site2 = barcode[23]
    site3 = barcode[24]
    site4 = barcode[25]
    site5 = barcode[26]
    site6 = barcode[27]
    site7 = barcode[28]
    site8 = barcode[29]
    site9 = barcode[30]
    site10 = barcode[31]
    if PASS != 'FAIL':
        output.write(str(cell_id) + '\t' + str(umi) + '\t' + str(PASS) + '\t' + str(site1) + '\t' + str(site2) + '\t' +str(site3) + '\t' +str(site4) + '\t' +str(site5) + '\t' +str(site6) + '\t' +str(site7) + '\t' +str(site8) + '\t' +str(site9) + '\t' +str(site10) + '\n')                     
output.close()
stats.close()

# Now let's make a barcode data frame that will merge well with Seurat objects.
print('Counting unique lineage barcodes...')
stats_brief = open(stats_file_summary)
stats_brief.readline()
cells_dict = {}

for line in stats_brief.readlines():
    line = line.split("\t")
    cell = line[0]
    PASS = line[2]
    barcode = str(line[3]+'_'+line[4]+'_'+line[5]+'_'+line[6]+'_'+line[7]+'_'+line[8]+'_'+line[9]+'_'+line[10]+'_'+line[11]+'_'+line[12])
    barcode = barcode.strip("\n")
    if cell in cells_dict: 
        cells_dict[cell].append(barcode)
    else: 
        cells_dict[cell] = [barcode]
        
stats_brief.close()


for cell in cells_dict.keys():
    unique_barcodes = set(cells_dict[cell])
    count = len(unique_barcodes)
    cells_dict[cell].append(count)


for cell in cells_dict: 
    test = cells_dict[cell]
    data = Counter(test)
    abundant = max(test, key=data.get)
    cells_dict[cell].append(abundant)

print('Writing "all barcodes" output file...')
output = open(output_all_barcodes, 'w')
output.write('cell.id' + '\t' 'barcode.count' + '\t' + 'barcodes' + '\t'+ 'site1'+ '\t'+ 
             'site2'+ '\t'+'site3'+ '\t'+'site4'+ '\t'+'site5'+ '\t'+'site6'+ '\t'+'site7'+ '\t'+'site8'+ '\t'+
             'site9'+ '\t'+'site10'+'\n')
for cell,barcode in cells_dict.items(): 
    if len(cells_dict[cell]) == 3: 
        sites = barcode[0].split("_")
        output.write(cell+'\t'+str(barcode[-2])+'\t'+barcode[0]+'\t'+sites[0]+'\t'+sites[1]+'\t'+sites[2]
                     +'\t'+sites[3]+'\t'+sites[4]+'\t'+sites[5]+'\t'+sites[6]+'\t'+sites[7]+'\t'+sites[8]
                     +'\t'+sites[9]+'\n')
    elif len(cells_dict[cell]) >=4:
        sites = barcode[0].split("_")
        output.write(cell+'\t'+str(barcode[-2])+'\t'+barcode[0]+'\t'+sites[0]+'\t'+sites[1]+'\t'+sites[2]
                     +'\t'+sites[3]+'\t'+sites[4]+'\t'+sites[5]+'\t'+sites[6]+'\t'+sites[7]+'\t'+sites[8]
                     +'\t'+sites[9]+'\n')
        for i in range(1,len(cells_dict[cell])-2): 
            sites = barcode[i].split("_")
            output.write(' '+'\t'+' '+'\t'+str(barcode[i])+'\t'+sites[0]+'\t'+sites[1]+'\t'+sites[2]
                     +'\t'+sites[3]+'\t'+sites[4]+'\t'+sites[5]+'\t'+sites[6]+'\t'+sites[7]+'\t'+sites[8]
                     +'\t'+sites[9]+'\n')
output.close()

print('Writing "top barcode" output file...')
output = open(output_top_barcodes, 'w')
output.write('cell.id' + '\t' 'barcode.count' + '\t' + 'barcodes' + '\t'+ 'site1'+ '\t'+ 
             'site2'+ '\t'+'site3'+ '\t'+'site4'+ '\t'+'site5'+ '\t'+'site6'+ '\t'+'site7'+ '\t'+'site8'+ '\t'+
             'site9'+ '\t'+'site10'+'\n')

for cell,barcode in cells_dict.items(): 
    sites = barcode[-1].split("_")
    output.write(cell+'\t'+str(barcode[-2])+'\t'+barcode[-1]+'\t'+sites[0]+'\t'+sites[1]+'\t'+sites[2]
                    +'\t'+sites[3]+'\t'+sites[4]+'\t'+sites[5]+'\t'+sites[6]+'\t'+sites[7]+'\t'+sites[8]
                     +'\t'+sites[9]+'\n')
output.close()

print("Complete!\n")
print(output_all_barcodes + " contains the complete list of barcodes associated with a cell id.\n")
print(output_top_barcodes + " contains the most abundant barcode associated with a cell id. \nUse this file for merging with your transcriptomic atlas data in R.\n")

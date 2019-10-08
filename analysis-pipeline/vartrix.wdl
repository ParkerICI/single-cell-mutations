workflow genomicsVarTrix {
  ## required files
  File cell_barcodes
  File bamFile
  File baiFile

  File vcfFile
  File vcfiFile

  ## reference files
  File refFasta
  File refFastai


  String sample_id


call VarTrix {
    input: 
      cell_barcodes=cell_barcodes,
      bamFile=bamFile,
      baiFile=baiFile,
      RefFasta=refFasta, 
      RefIndex=refFastai,
      vcfFile=vcfFile,
      vcfiFile=vcfiFile,
      sample_name=sample_id
  }
}

task VarTrix {

  File cell_barcodes
  File bamFile
  File baiFile
  File RefFasta
  File RefIndex
  File vcfFile
  File vcfiFile
  String sample_name


  command <<<
    vartrix -v ${vcfFile} 
              -b ${bamFile} 
              -f ${RefFasta} 
              -c ${cell_barcodes}
              -o ${sample_name}

    vawk '{print $1,$2}' vcfFile > SNV.loci.txt
    sed -i 's/\s/:/g' SNV.loci.txt 

  >>>

  output {
    File output_file = "${sample_name}"
    File snv_loci = "SNV.loci.txt"
  }

}





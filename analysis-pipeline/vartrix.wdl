workflow genomicsVarTrix {
  ## required files
  File cell.barcodes
  File bamFile
  File baiFile

  File vcfFile
  File vcfiFile

  ## reference files
  File refFasta
  File refFastai


  String output.path=sample.id


call VarTrix {
    input: 
      cell.barcodes=cell.barcodes,
      bamFile=bamFile,
      baiFile=baiFile,
      RefFasta=refFasta, 
      RefIndex=refFastai,
      vcfFile=vcfFile,
      vcfFilei=vcfFilei,
      sample_name=output.path
  }
}

task VarTrix {

  File cell.barcodes
  File bamFile
  File bamIndex
  File RefFasta
  File RefIndex
  File vcfFile
  File vcfFilei
  String sample_name


  command {
    vartrix -v ${vcfFile} 
              -b ${bamFile} 
              -f ${refFasta} 
              -c ${cell.barcodes}
              -o ${sample_name}

    vawk '{print $1,$2}' vcfFile > SNV.loci.txt
    sed -i 's/\s/:/g' SNV.loci.txt 
  }

  output {
    File output = "${sample_name}.txt"
    File snv.loci = SNV.loci.txt
  }

}





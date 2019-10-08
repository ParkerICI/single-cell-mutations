workflow 10x.genomicsVarTrix {
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

  # Runtime
  String docker

call GetBwaVersion {
    input:
      docker = docker
  }


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

  File RefFasta
  File RefIndex
  File RefDict
  String sample_name
  File bamFile
  File bamIndex

  command {
    vartrix -v ${vcfFile} 
              -b ${bamFile} 
              -f ${refFasta} 
              -c ${cell.barcodes}
              -o ${sample_name}
  }

  output {
    File output = "${sample_name}.txt"
  }

  command { 
     vawk '{print $1,$2}' vcfFile > SNV.loci.txt
     sed -i 's/\s/:/g' SNV.loci.txt 
   }

  ## this needs toe be fixed.
   output { 
     File snv.loci = SNV.loci.txt} 




}





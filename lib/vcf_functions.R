vcf_to_tidy<-function(input_vcfR,happy_individual=NA){
  message(glue('Converting VCF to tidy form..'))
  vcf_tidy<-vcfR2tidy(input_vcfR)
  gt_df<-vcf_tidy$gt
  if (!is.na(happy_individual)){
    gt_df<-gt_df%>%filter(Indiv==happy_individual)
  }
  info_df<-vcf_tidy$fix
  if (nrow(info_df)!=nrow(gt_df)){stop(glue('The number of rows in the fix part of the vcf is different from the number in the gt part'))}
  vcf_df<-info_df%>%
    bind_cols(gt_df%>%select(-c(ChromKey,POS,Indiv)))
  return(vcf_df)
}

gt_merge<-gt_df%>%filter(Indiv=='QUERY')%>%select(gt_BD,gt_BLT,gt_BVT)%>%
  bind_cols(gt_df%>%filter(!Indiv%in%c('QUERY','TRUTH'))%>%select(-c(gt_BD,gt_BLT,gt_BVT)))%>%
  bind_cols(vcf_tidy$fix%>%select(-c(ChromKey,POS)))
gt_merge<-gt_merge%>%filter(gt_BLT!='nocall')
z<-gt_merge%>%filter(is.na(gt_GT))

filter_happy_vcf_by_region<-function(happy_vcf_file_name,regions_of_interest){
  open_file_command <-ifelse(grepl('.gz',happy_vcf_file_name),'gunzip -c','cat')
  regions_query<-paste0(regions_of_interest,collapse='\\|')
  happy_file_name_base<-basename(happy_vcf_file_name)%>%str_replace_all('.(vcf|gz)','')
  happy_file_folder<-happy_vcf_file_name%>%str_replace('/[^/]+$','')
  out_file_name<-glue('{happy_file_name_base}.regions_of_interest.vcf')
  filter_command<-glue('{open_file_command} {happy_vcf_file_name} | grep "#\\|{regions_query}" > {happy_file_folder}/{out_file_name}')
  message(glue('Running {filter_command}..'))
  system(filter_command)
  return(as.character(glue('{happy_file_folder}/{out_file_name}')))
}

normalize_vcf<-function(vcf_file_name,reference_file,output_folder){
  norm_file_name<-stringr::str_replace_all(basename(vcf_file_name),'.(vcf|gz)','')
  norm_file_name<-glue('{output_folder}/{norm_file_name}.norm.vcf.gz')
  normalization_command<-glue('bcftools norm --do-not-normalize {vcf_file_name} -O z -m - -f {reference_file} > {norm_file_name}')
  message(glue('Normalizing VCF: {normalization_command}'))
  system(normalization_command)
  return(norm_file_name)
}

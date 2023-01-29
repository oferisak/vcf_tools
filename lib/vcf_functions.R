vcf_to_tidy<-function(input_vcfR,happy_individual=NA){
  vcf_tidy<-vcfR2tidy(input_vcfR)
  gt_df<-vcf_tidy$gt
  info_df<-vcf_tidy$fix
  vcf_df<-gt_df%>%left_join(info_df)
  if (!is.na(happy_individual)){
    vcf_df<-vcf_df%>%filter(Indiv==happy_individual)
  }
  return(vcf_df)
}

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

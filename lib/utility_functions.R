# a function that takes as input a string with a condition (e.g "QUAL>20") and returns it as a condition
get_logical_vector <- function(df, condition) {
  condition_to_parse<-paste0(paste("df$",condition,sep = ''),collapse='&')
  print(condition_to_parse)
  logical_vector <- eval(parse(text = condition_to_parse))
  return(logical_vector)
}

# a function that takes in a list of conditions and annotates the given input table according to the given conditions
bin_preds<-function(input_table,conditions,col_name='bin_pred',default='intermediate'){
  input_table[,col_name]<-default
  for (i in 1:length(conditions)){
    condition_name<-names(conditions[i])
    condition=conditions[[i]]
    matching_rows<-get_logical_vector(input_table, condition)
    input_table[matching_rows,col_name]<-condition_name
  }
  return(input_table)
}

calculate_bin_counts<-function(test_data,conditions,preds_col,by_vartype=F,output_ppv=F){
  test_data_fixed<-bin_preds(test_data,col_name=preds_col,conditions)
  if (by_vartype){
    count_data_tree<-test_data_fixed%>%count(!!sym(preds_col),gt_BD,gt_BVT)
    totals<-count_data_tree%>%group_by(!!sym(preds_col),gt_BVT)%>%summarize(total=sum(n))
  }else{
    count_data_tree<-test_data_fixed%>%count(!!sym(preds_col),gt_BD)
    totals<-count_data_tree%>%group_by(!!sym(preds_col))%>%summarize(total=sum(n))
  }
  count_data_tree<-count_data_tree%>%
    left_join(totals)%>%
    mutate(rate=n/total)
  
  if (output_ppv & by_vartype){
    count_data_tree<-
      count_data_tree%>%select(-rate)%>%
      pivot_wider(names_from = gt_BD,values_from = n)%>%
      left_join(count_data_tree%>%group_by(gt_BD,gt_BVT)%>%summarize(total_in_BD=sum(n))%>%filter(gt_BD=='TP'))%>%
      group_by(!!sym(preds_col),gt_BVT)%>%
      summarize(PPV=TP/(TP+FP),
                sensitivity=TP/total_in_BD)
  }
  if (output_ppv & !by_vartype){
    count_data_tree<-
      count_data_tree%>%select(-rate)%>%
      pivot_wider(names_from = gt_BD,values_from = n)%>%
      bind_cols(count_data_tree%>%group_by(gt_BD)%>%summarize(total_in_BD=sum(n))%>%filter(gt_BD=='TP'))%>%
      group_by(!!sym(preds_col))%>%
      summarize(PPV=TP/(TP+FP),
                sensitivity=TP/total_in_BD)
  }
  
  return(count_data_tree)
  
}


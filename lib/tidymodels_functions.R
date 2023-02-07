# Split data into train and test
split_train_test<-function(ml_data,strata,seed=123,train_prop=0.75){
  set.seed(seed)
  split <- ml_data%>%rsample::initial_split(strata = strata,prop = train_prop)
  train_test<-list()
  train_test[['train']] <- rsample::training(split)
  train_test[['test']] <- rsample::testing(split)
  return(train_test)
}

# Setup what model parameters you want to tune during training
setup_decision_tree_tuning<-function(engine='rpart'){
  tune_tree <- 
    decision_tree(
      cost_complexity = tune(), #tune() is a placeholder for an empty grid 
      tree_depth = tune() #we will fill these in the next section
    ) %>% 
    set_engine("rpart") %>% 
    set_mode("classification") 
  return(tune_tree)  
}

# Workflow setup
setup_decision_tree_workflow<-function(recipe,tuning_setup){
  decision_tree_wf <- workflow() %>% 
    add_model(tuning_setup) %>% 
    add_recipe(recipe)
  return(decision_tree_wf)
}

# Hyperparameter evaluation
train_model_on_tuning_grid<-function(workflow_with_tune,vfold,selected_grid,
                                     metrics=yardstick::metric_set(accuracy, roc_auc,mn_log_loss)){
  tuning_result <- 
    workflow_with_tune %>% 
    tune::tune_grid(
      resamples = vfold, #This is the 10 fold cross validation variable we created earlier
      grid = selected_grid, #This is the tuning grid
      metrics = metrics
    )
  return(tuning_result)
}

train_decision_tree<-function(rec,vfold,tuning_grid,metric_for_model_selection='roc_auc'){
  # set up the tuning process
  decision_tree_tuning_setup <- setup_decision_tree_tuning()
  # set up the workflow
  decision_tree_wf <- setup_decision_tree_workflow(recipe = rec,
                                                   tuning_setup = decision_tree_tuning_setup)
  
  # Evaluate model performance on cross validation using the different tuning parameters
  tuning_res<-train_model_on_tuning_grid(workflow_with_tune = decision_tree_wf,
                                         vfold = vfold,
                                         selected_grid = tuning_grid)
  
  # collected metrics
  collected_mets <- tune::collect_metrics(tuning_res)
  print(collected_mets)
  
  best_tree <- tuning_res %>% 
    tune::select_best(metric_for_model_selection)
  
  final_decision_tree_wf <- 
    decision_tree_wf %>% 
    finalize_workflow(best_tree) #Finalise workflow passes in our best tree
  
  print(final_decision_tree_wf)
  
  final_decision_tree <- 
    final_decision_tree_wf %>% 
    fit(data = train_data)
  
  print(final_decision_tree)
  return(final_decision_tree)
}

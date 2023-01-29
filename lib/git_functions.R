library(usethis)

setup_git<-function(){
  ## now, working inside "the package", initialize git repository
  use_git()
}

create_git_token<-function(){
  ## create github token
  create_github_token()
  # set it up
  gitcreds::gitcreds_set()
}

create_git_repo<-function(){
  ## create github repository and configure as git remote
  use_github()
}
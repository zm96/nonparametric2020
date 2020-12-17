VaryCurve=function(u,mod=2)
{
  if(mod==2)
  {
    res=rbind(u^2/2,
              exp((-(u>0)*1.44-(u<0)/1.44)*u^2/4),
              exp(-0.72^(-1)*(u+1)^2)+0.7*exp(-0.98^(-1)*(u-1)^2),
              (0.4*sin(pi*(u-0.5)/2.5)+0.8)/exp((u+(u+0.5)^2*I(u>-0.5))/6))
  }
  if(mod==1)
  {
    res=rbind(u^3/8,
              sin(pi*u/2),
              u*(1-u)/4,
              exp(-u^2))
  }
  return(res)
}

# 
# VCM=function(u_grid,Y,U,X,W_sample,W_cv,start_ind,end_ind,kernel_type=2)
# {
#   CV.object<-function(h) #objective function of cross-validation selection of h
#   {
#     res=VCM_cv_cpp(X,U,Y,h,W_sample,W_cv,start_ind,end_ind,kernel_type)
#     return(res)
#   }
#   h=optimize(CV.object,interval=c(0.01,1))$minimum 
#   fit_ind=VCM_cpp(X,U,Y,h,u_grid,W_sample,kernel_type)
#   return(list(h=h,fit=fit_ind[1:(ncol(X)),]))
# }

VCM=function(u_grid,Y,U,X,W,h_range,kernel_type=2)
{
  W_cv=W
  W_sample=W
  start_ind=seq(1,length(Y)) #starting index of each pool
  end_ind=seq(1,length(Y))   #ending index of each pool
  CV.object<-function(h) #objective function of cross-validation selection of h
  {
    res=VCM_cv_cpp(X,U,Y,h,W_sample,W_cv,start_ind,end_ind,kernel_type)
    return(res)
  }
  h=optimize(CV.object,interval=h_range)$minimum 
  fit=VCM_cpp(X,U,Y,h,u_grid,W_sample,kernel_type)
  return(list(h=h,fit=fit[1:ncol(X),]))
}
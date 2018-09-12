lexEvalTypeToString <- function(func_index, eval_funcs) {
  #cosine
  if (func_index - length(eval_funcs[[1]]) <= 0) {
    return("$K_{Cos}$")
  }
  
  func_index <- func_index - length(eval_funcs[[1]])
  
  #polynomial_kernel_func
  if (func_index - length(eval_funcs[[2]]) <= 0) {
    return("$K_{Poly}$")
  }
  
  func_index <- func_index - length(eval_funcs[[2]])
  
  #p_spectrum_kernel_func
  if (func_index - length(eval_funcs[[3]]) <= 0) {
    return("$K_{Spec}$")
  }
  
  func_index <- func_index - length(eval_funcs[[3]])
  
  #constant_kernel_func
  if (func_index - length(eval_funcs[[4]]) <= 0) {
    return("$K_{Const}$")
  }
  
  return("NOT FOUND")
  
}

cfgEvalTypeToString <- function(func_index, eval_funcs) {
  #  cfgsim_kernel_funcs <- unlist(list(exponential_diffusion_kernel_func, laplacian_exponential_diffusion_kernel_func))

    #exponential_diffusion_kernel_func
  if (func_index - length(eval_funcs[[1]]) <= 0) {
    return("$K_{Exp}$")
  }
  
  func_index <- func_index - length(eval_funcs[[1]])
  
  #laplacian_exponential_diffusion_kernel_func
  if (func_index - length(eval_funcs[[2]]) <= 0) {
    return("$K_{LED}$")
  }
 
  return("NOT FOUND")
  
}

freqEvalTypeToString <- function(func_index, eval_funcs) {
  #    freqsim_kernel_funcs <- unlist(list(polynomial_kernel_func, gaussian_kernel_func))
  
  #polynomial_kernel_func
  if (func_index - length(eval_funcs[[1]]) <= 0) {
    return("$K_{Poly}$")
  }
  
  func_index <- func_index - length(eval_funcs[[1]])
  
  #gaussian_kernel_func
  if (func_index - length(eval_funcs[[2]]) <= 0) {
    return("$K_{Gaus}$")
  }
  
  return("NOT FOUND")
  
}

#input parameters are unlisted so indices of all-funcs can be mapped to parameters
evalFuncToString = function(func_index, type, parameters) {
  
  eval_string <-  switch (type,
                          "cfg" = cfgEvalTypeToString(func_index, parameters),
                          "freq" = freqEvalTypeToString(func_index, parameters),
                          "lex" = lexEvalTypeToString(func_index, parameters)
  )
  
  parameters <- unlist(parameters)
  par_string = parameters[func_index]
  
  if (par_string != 0){
    eval_string <- paste(eval_string, " (", "$", par_string, "$)", sep="")
  }
  
  return (eval_string)
}
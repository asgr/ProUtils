quan_wt_array = function(array, probs = 0.5, wt = NULL, type = 7){
  if(is.list(array)){
    dim_array = c(dim(array[[1]]), length(array))
    array = unlist(array)
  }else{
    dim_array = dim(array)
  }

  dim(array) = c(dim_array[1] * dim_array[2], dim_array[3])

  if(!is.null(wt)){
    if(is.list(wt)){
      dim_wt = c(dim(wt[[1]]), length(wt))
      wt = unlist(wt)
    }else{
      dim_wt = dim(wt)
    }

    if(!all(dim_array == dim_wt)){
      stop('Dimensions of array and wt do not match!')
    }

    dim(wt) = c(dim_wt[1] * dim_wt[2], dim_wt[3])
  }

  output = quan_wt_mat_row(array, probs, wt, type)

  if(length(probs) == 1){
    dim(output) = dim_array[1:2]
  }else{
    dim(output) = c(dim_array[1:2], length(probs))
  }

  return(output)
}

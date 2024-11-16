#' @title slot accessor functions for class 'tprtrees'
#' @description get and set slot values
#' @param obj object of class 'tprtrees'
#' @details Getting and setting the values of the different slots of
#' 'tprTrees'-objects. For slot \code{mono} no setting function has been defined,
#' as this slot is computed by \code{\link{check_monotonicity}} and should not
#' be reset by users.
#' @return the accessor functions return the value of the specified slot and
#' the setting functions update the object
#' @name Slot_accessors
#' @seealso \code{\link{tprTrees-class}}, \code{\link{tprTrees}}
NULL

#' @describeIn Slot_accessors getting slot 'spp' of \code{obj}
#' @export
setGeneric("spp",
           function(obj) standardGeneric("spp"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("spp", signature = "tprTrees", function(obj){
            return(obj@spp)
          })

#' @describeIn Slot_accessors setting 'spp' slot of object
#' @details Setting of \code{spp} requires mode integer. For convenience,
#' \code{value} is coerced by \code{\link{as.integer}}.
#'
#' Setting \code{spp} and \code{H}, a vector of length equal
#' \code{length(spp(obj))} is required.
#' @export
setGeneric("spp<-", function(obj, value) standardGeneric("spp<-"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("spp<-", signature = "tprTrees", function(obj, value) {
  if(length(value) < length(obj)){
    value <- rep(value[1], length(obj))
  }
  obj@spp <- as.integer(value)
  obj
})

#' @describeIn Slot_accessors getting slot 'Dm' of \code{obj}
#' @export
setGeneric("Dm",
           function(obj) standardGeneric("Dm"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("Dm", signature = "tprTrees",
          function(obj){
            x <- obj@Dm
            if(length(obj@spp) == 1){
              return(x[[1]])
            } else {
              return(x)
            }
          })

#' @param value depending on slot, see details
#' @describeIn Slot_accessors setting 'Dm' slot of object
#' @details For setting slots \code{Dm} and \code{Hm} \code{value} must be a
#' list of vectors of length equal \code{length(spp(obj))} and the length of
#' each vector must correspond to the length of the vectors in \code{Hm} and
#' \code{Dm}.
#' @export
setGeneric("Dm<-", function(obj, value) standardGeneric("Dm<-"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("Dm<-", signature = "tprTrees", function(obj, value) {
  obj@Dm <- value
  obj
})

#' @describeIn Slot_accessors getting slot 'Hm' of \code{obj}
#' @export
setGeneric("Hm",
           function(obj) standardGeneric("Hm"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("Hm", signature = "tprTrees",
          function(obj){
            x <- obj@Hm
            if(length(obj@spp) == 1){
              return(x[[1]])
            } else {
              return(x)
            }
          })

#' @param value depending on slot, see details
#' @describeIn Slot_accessors setting 'Hm' slot of object
#' @export
setGeneric("Hm<-", function(obj, value) standardGeneric("Hm<-"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("Hm<-", signature = "tprTrees", function(obj, value) {
  obj@Hm <- value
  obj
})

#' @describeIn Slot_accessors getting slot 'Ht' of \code{obj}
#' @export
setGeneric("Ht", function(obj) standardGeneric("Ht"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("Ht", signature = "tprTrees",
          function(obj){
            return(obj@Ht)
          })

#' @describeIn Slot_accessors setting 'Ht' slot of object
#' @export
setGeneric("Ht<-", function(obj, value) standardGeneric("Ht<-"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("Ht<-", signature = "tprTrees", function(obj, value) {
  if(length(value) < length(obj)){
    value <- rep(value[1], length(obj))
  }
  obj@Ht <- value
  obj
})

#' @describeIn Slot_accessors getting slot 'sHt' of \code{obj}
#' @export
setGeneric("sHt", function(obj) standardGeneric("sHt"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("sHt", signature = "tprTrees",
          function(obj){
            return(obj@sHt)
          })

#' @describeIn Slot_accessors setting 'sHt' slot of object
#' @export
setGeneric("sHt<-", function(obj, value) standardGeneric("sHt<-"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("sHt<-", signature = "tprTrees", function(obj, value) {
  if(length(value) < length(obj)){
    value <- rep(value[1], length(obj))
  }
  obj@sHt <- value
  obj
})

#' @describeIn Slot_accessors getting slot 'monotone' of \code{obj}
#' @export
setGeneric("mono",
           function(obj) standardGeneric("mono"))

#' @describeIn Slot_accessors method for class 'tprTrees'
setMethod("mono", signature = "tprTrees",
          function(obj){
            return(obj@monotone)
          })

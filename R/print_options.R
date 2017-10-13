#' List all control options in the \code{cvma} package.
#' 
#' @export
#' @examples
#' list_control_options()
list_control_options <- function(){
    message("All ensemble functions in cvma:\n")
    print(ls("package:cvma", pattern = "ensemble_"))
	
	message("All weight functions used to select weighted super learner in cvma:\n")
    print(ls("package:cvma", pattern = "weight_sl"))	

    message("All weight functions used to select weighted outcomes in cvma:\n")
    print(ls("package:cvma", pattern = "weight_y"))

	message("All risk functions used to select weighted super learner in cvma:\n")
    print(ls("package:cvma", pattern = "optim_risk_sl"))	

    message("All risk functions used to select weighted outcomes in cvma:\n")
    print(ls("package:cvma", pattern = "optim_risk_y"))

	message("All cv-risk functions used to evaluate super learner predictions of single outcomes in cvma:\n")
    print(ls("package:cvma", pattern = "cv_risk_sl"))	

    message("All cv-risk functions used to evaluate super learner predictions of combined outcomes in cvma:\n")
    print(ls("package:cvma", pattern = "cv_risk_y"))
}
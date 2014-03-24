zhome <- "$FOLDER$"

run <- function(files_to_run)
{
	zregion   <- "$REGION$"
	folder    <- "RUX2_region$SUBREGION$"
	filetorun <- c("continue.r", "mcmc.r","analysis.r","outbreaks_text.r","outbreaks_pdf.r","outbreaks_temporal.r","outbreaks_kml.r","outbreaks_map.r")

	for (i in files_to_run)
	{
		source_file <- paste(zhome, zregion, "/", folder, "/", filetorun[i], sep="")
		cat("\n\n========================================\nrunning ",source_file,"\n=========================================\n\n")
		source(source_file)
	}
}

continue <- function()
{
  run(c(1,3,5,7,8))
}

analysis <- function()
{
  run(c(3,5,7,8))
}

full_model <- function()
{
  run(c(2,3,5,7,8))
}

cat("Welcome to EpiclustR\n")
cat("====================\n\n")
cat("Model run from", $BEGINTIME$, "to", $ENDTIME$,"\n\n")
cat("This is the simple R interface for re-running portions of the model.\n\n")
cat("Type continue() to continue a model that didn't finish all iterations\n")
cat("Type analysis() to run the analysis stage of the model\n")
cat("Type full_model() to run the full model from scratch (overriding any progress)\n")


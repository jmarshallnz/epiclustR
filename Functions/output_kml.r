
#
# Write a KML file header
#

write_kml_header <- function(file, name)
{
  cat(file=file, sep="", append=F, "<?xml version=\"1.0\" encoding=\"us-ascii\"?>\n")
  cat(file=file, sep="", append=T, "<kml xmlns=\"http://earth.google.com/kml/2.1\">\n")
  cat(file=file, sep="", append=T, " <Document>\n")
  cat(file=file, sep="", append=T, "  <name>", name, "</name>\n")
}

write_kml_footer <- function(file)
{
  cat(file=file, sep="", append=T, " </Document>\n</kml>\n")
}

week_year <- function(t)
{
  return(sprintf("%.0f, %.0f", floor(as.POSIXlt(t)$yday / 7) + 1, as.POSIXlt(t)$year + 1900));
}

write_epidemic_header <- function(file, begin, region, num_cases, prob)
{
  description <- sprintf("Week: %s<br>Beginning: %s<br>Number of cases: %i<br>Probability: %f", week_year(begin), begin, num_cases, prob);
  name <- sprintf("Week %s, %s", week_year(begin), region);
  cat(file=file, sep="", append=T, "  <Folder>\n   <name>", name, "</name>\n")
  cat(file=file, sep="", append=T, "   <description><![CDATA[", description, "]]></description>\n")
}

write_epidemic_footer <- function(file)
{
  cat(file=file, sep="", append=T, "  </Folder>\n")
}

write_epidemic_case <- function(file, begin, end, long, lat, case_id, mb_id, prob, color, num_cases)
{
  description <- sprintf("Week: %s<br> Beginning: %s<br>Longitude: %f</br>Latitude: %f<br>Cases: %s<br>Spatial Unit: %s<br>Probability: %f", week_year(begin), begin, long, lat, case_id, mb_id, prob);
  scale <- sqrt(num_cases);

  cat(file=file, sep="", append=T, "   <Placemark>\n")
  cat(file=file, sep="", append=T, "    <description><![CDATA[", description, "]]></description>\n")

  cat(file=file, sep="", append=T, "    <TimeSpan>\n")
  cat(file=file, sep="", append=T, "     <begin>", as.character(begin), "</begin>\n")
  cat(file=file, sep="", append=T, "     <end>", as.character(end), "</end>\n")
  cat(file=file, sep="", append=T, "    </TimeSpan>\n")

  cat(file=file, sep="", append=T, "    <Point>\n")
  cat(file=file, sep="", append=T, "     <coordinates>", long, ",", lat, ",0</coordinates>\n")
  cat(file=file, sep="", append=T, "    </Point>\n")

  cat(file=file, sep="", append=T, "    <Style>\n")
  cat(file=file, sep="", append=T, "     <IconStyle>\n")
  cat(file=file, sep="", append=T, "      <Icon>http://maps.google.com/mapfiles/kml/shapes/shaded_dot.png</Icon>\n")
  cat(file=file, sep="", append=T, "      <color>", color, "</color>\n")
  cat(file=file, sep="", append=T, "      <scale>", scale, "</scale>\n")
  cat(file=file, sep="", append=T, "     </IconStyle>\n")
  cat(file=file, sep="", append=T, "    </Style>\n")

  cat(file=file, sep="", append=T, "   </Placemark>\n")
}

#write_kml_header("test.kml", "test output")

#write_spatial_header("test.kml", "region1")
#write_epidemic("test.kml", "2001-04-20", "2001-05-19", 170.65, -47.36, 5, 0.5);
#write_epidemic("test.kml", "2001-06-20", "2001-07-19", 170.65, -47.36, 5, 0.5);
#write_epidemic("test.kml", "2001-08-20", "2001-09-19", 170.65, -47.36, 5, 0.5);
#write_spatial_footer("test.kml")

#write_spatial_header("test.kml", "region1")
#write_epidemic("test.kml", "2001-05-20", "2001-06-19", 170.85, -47.36, 5, 0.5);
#write_epidemic("test.kml", "2001-07-20", "2001-08-19", 170.85, -47.36, 5, 0.5);
#write_epidemic("test.kml", "2001-09-20", "2001-10-19", 170.85, -47.36, 5, 0.5);
#write_spatial_footer("test.kml")

#write_kml_footer("test.kml")

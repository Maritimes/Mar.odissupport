require(RODBC)
channel = make_oracle_cxn()[[2]]
require(arcgisbinding)
arc.check_product()

serviceMaker(channel, write.path ="C:/Users/mcmahonm/Documents/ArcGIS/rPED.gdb")

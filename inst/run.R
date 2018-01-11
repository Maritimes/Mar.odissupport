channel = Mar.utils::make_oracle_cxn()[[2]]
arcgisbinding::arc.check_product()

serviceMaker(channel, write.path ="C:/Users/mcmahonm/Documents/ArcGIS/rPED.gdb")

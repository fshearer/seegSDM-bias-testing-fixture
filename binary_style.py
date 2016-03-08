import mapnik
import sys

diseases = ("cchf", "chik", "deng", "hat", "melio", "nwcl", "owcl", "scrub")
variants = ("1A", "1B", "1C", "2C", "2D", "3B", "4A", "5A", "6A")

for disease in diseases:
  for variant in variants:
    geotiff = "all/results/"+disease+"_"+variant+"/binary_prediction.tif"
    png = "all/results/"+disease+"_"+variant+"/binary_prediction.png"
    symbolizer = mapnik.RasterSymbolizer()
    colorizer = mapnik.RasterColorizer(default_mode=mapnik.COLORIZER_LINEAR, default_color=mapnik.Color("transparent"))
    symbolizer.colorizer = colorizer
    colorizer.add_stop(-9999.0, mapnik.Color("#000000"))
    colorizer.add_stop(0.0, mapnik.Color("#eaeaea"))
    colorizer.add_stop(1.0, mapnik.Color("#454545"))
    style = mapnik.Style()
    rule = mapnik.Rule()
    rule.symbols.append(symbolizer)
    style.rules.append(rule)
    layer = mapnik.Layer("raster")
    layer.datasource = mapnik.Gdal(file=geotiff, nodata=-9999, band=1)
    layer.styles.append("raster")
    world = mapnik.Map(1656, 667)
    world.append_style("raster", style)
    world.layers.append(layer)
    world.zoom_all()
    mapnik.render_to_file(world, png, 'png')

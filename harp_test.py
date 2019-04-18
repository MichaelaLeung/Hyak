import harp
# Create a product in Python and export it as a NetCDF file.
product = harp.Product()

# Import an S5P L2 HCHO product.
hcho_product = harp.import_product("SCI_NL__1PYDPA20020803_003946_000057692008_00160_02218_0000.N1",
                                   "solar_zenith_angle < 60 [degree]; latitude > 30 [degree_north]; latitude < 60 [degree_north]")

# Pretty print information about the product.
print(hcho_product)



# Create a product in Python and export it as a NetCDF file.
product = harp.Product()
harp.export_product(product, "SCI_NL__1PYDPA20120408_002516_000060013113_00261_52861_0000.N1")

# Add some variables to the product.
product.foo = harp.Variable("foo")
product.strings = harp.Variable(numpy.array(("foo", "bar", "baz")), ["time"])
product.temperature = harp.Variable(numpy.ones((3, 5), dtype=numpy.float32),
                                    ["time", None])
product.temperature.unit = "K"
product.temperature.description = "temperature"

# Pretty print information about the product.
print(product)

# Pretty print information about the variable 'temperature'.
print(product.temperature)

# Set valid minimum value of the variable 'temperature'. Note the use of item
# access syntax instead of attribute access syntax.
product["temperature"].valid_min = 0.0
print(product.temperature)

# Export the updated product as an HDF4 file.
harp.export_product(product, "non-empty.hdf", file_format="hdf4")

# Convert the product to an OrderedDict.
dict_product = harp.to_dict(product)

# Import an S5P L2 HCHO product.
hcho_product = harp.import_product("S5P_NRTI_L2__HCHO___20080808T224727_20080808T234211_21635_01_021797_00000000T000000.nc",
                                   "solar_zenith_angle < 60 [degree]; latitude > 30 [degree_north]; latitude < 60 [degree_north]")

# Pretty print information about the product.
print(hcho_product)

# Export the product as a HARP compliant data product.
harp.export_product(hcho_product, "hcho.h5", file_format='hdf5', hdf5_compression=6)


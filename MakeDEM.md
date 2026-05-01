## DEMデータ変換の例
### 【日本】
250mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 123.0000 46.000 146.0000 24.0000 ^
  -tr 0.0027728 0.0022533 ^
  -co COMPRESS=LZW -co TILED=YES ^
  U:\DEM\all_japan.vrt U:\DEM\Japan_250m.tif
```
50mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 123.0000 46.000 146.0000 24.0000 ^
  -tr 0.00055456 0.00045066 ^
  -co COMPRESS=LZW -co TILED=YES ^
  U:\DEM\all_japan.vrt U:\DEM\Japan_50m.tif
```
### 【北海道】
250mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 139.0000 46.000 146.0000 41.3333 ^
  -tr 0.0031707 0.0022495 ^
  -co COMPRESS=LZW -co TILED=YES ^
  U:\DEM\all_japan.vrt U:\DEM\hokkaidou_250m.tif
```
100mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 139.0000 46.000 146.0000 41.3333 ^
  -tr 0.001268 0.000900 ^
  -co COMPRESS=LZW -co TILED=YES ^
  U:\DEM\all_japan.vrt U:\DEM\hokkaidou_100m.tif
```
50mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 139.0000 46.000 146.0000 41.3333 ^
  -tr 0.0006341 0.0004499 ^
  -co COMPRESS=LZW -co TILED=YES ^
  U:\DEM\all_japan.vrt U:\DEM\hokkaidou_50m.tif
```
10mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 139.0000 46.000 146.0000 41.3333 ^
  -tr 0.0001268 0.0000900 ^
  -co COMPRESS=LZW -co TILED=YES ^
  U:\DEM\all_japan.vrt U:\DEM\hokkaidou_10m.tif
```

### 【長野】
250mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 36.70100 138.50000 35.80000 ^
  -tr 0.0027728 0.0022533 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Nagano_250m.tif
```
100mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 36.70100 138.50000 35.80000 ^
  -tr 0.0011092 0.0009014 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Nagano_100m.tif
```
50mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 36.70100 138.50000 35.80000 ^
  -tr 0.0005546 0.0004507 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Nagano_50m.tif
```
25mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 36.70100 138.50000 35.80000 ^
  -tr 0.00027728 0.00022533 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Nagano_25m.tif
```
10mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 36.70100 138.50000 35.80000 ^
  -tr 0.0001109 0.0000901 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Nagano_10m.tif
```

### 【群馬】
250mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 138.50000 37.10100 139.60900 36.20000 ^
  -tr 0.0027728 0.0022533 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma_250m.tif
```
100mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 138.50000 37.10100 139.60900 36.20000 ^
  -tr 0.0011092 0.0009014 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma_100m.tif
```
50mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 138.50000 37.10100 139.60900 36.20000 ^
  -tr 0.0005546 0.0004507 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma_50m.tif
```
25mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 138.50000 37.10100 139.60900 36.20000 ^
  -tr 0.00027728 0.00022533 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma_25m.tif
```
10mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 138.50000 37.10100 139.60900 36.20000 ^
  -tr 0.0001109 0.0000901 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma_10m.tif
```

### 【群馬～長野～富山】
250mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 37.10100 139.60900 35.80000 ^
  -tr 0.0027728 0.0022533 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma-Toyama_250m.tif
```
100mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 37.10100 139.60900 35.80000 ^
  -tr 0.0011092 0.0009014 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma-Toyama_100m.tif
```
50mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 37.10100 139.60900 35.80000 ^
  -tr 0.0005546 0.0004507 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma-Toyama_50m.tif
```
25mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 137.39100 37.10100 139.60900 35.80000 ^
  -tr 0.00027728 0.00022533 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Gunma-Toyama_25m.tif
```

### 【九州】
250mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 128.5430 34.015 132.1250 30.950 ^
  -tr 0.0026731 0.0022533 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Kyushu_250m.tif
```
100mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 128.5430 34.015 132.1250 30.950 ^
  -tr 0.00106924 0.00090132 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Kyushu_100m.tif
```
50mメッシュ
```
gdal_translate -of GTiff -r average ^
  -projwin 128.5430 34.015 132.1250 30.950 ^
  -tr 0.00053462 0.00045066 ^
  -co COMPRESS=LZW ^
  U:\DEM\all_japan.vrt U:\DEM\Kyushu_50m.tif
```

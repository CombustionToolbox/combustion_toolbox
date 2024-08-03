mkdir icon.iconset
sips -z 16 16       logo_CT_noversion.png --out icon.iconset/icon_16x16.png
sips -z 32 32       logo_CT_noversion.png --out icon.iconset/icon_16x16@2x.png
sips -z 32 32       logo_CT_noversion.png --out icon.iconset/icon_32x32.png
sips -z 64 64       logo_CT_noversion.png --out icon.iconset/icon_32x32@2x.png
sips -z 128 128     logo_CT_noversion.png --out icon.iconset/icon_128x128.png
sips -z 256 256     logo_CT_noversion.png --out icon.iconset/icon_128x128@2x.png
sips -z 256 256     logo_CT_noversion.png --out icon.iconset/icon_256x256.png
sips -z 512 512     logo_CT_noversion.png --out icon.iconset/icon_256x256@2x.png
sips -z 512 512     logo_CT_noversion.png --out icon.iconset/icon_512x512.png
sips -z 1024 1024   logo_CT_noversion.png --out icon.iconset/icon_512x512@x2.png
iconutil -c icns icon.iconset
rm -R icon.iconset

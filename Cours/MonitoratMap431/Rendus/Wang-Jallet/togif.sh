magick mogrify -format png -background white -alpha remove *.eps
magick *.png -delay 10 -loop 0 *.png deformation.gif
rm *.png

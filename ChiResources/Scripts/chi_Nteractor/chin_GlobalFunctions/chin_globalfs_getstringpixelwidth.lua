--================================== Get String Pixel Width
--Obtains the pixel width of a string of characters.
function chinGetStringPixelWidth(str)
    pixelWidth=0;
    for k=1,(string.len(str)) do
        charCode=string.byte(string.sub(str,k,k));
        pixelWidth=pixelWidth+chiGetCharacterWidth(charCode);
    end
    return pixelWidth;
end
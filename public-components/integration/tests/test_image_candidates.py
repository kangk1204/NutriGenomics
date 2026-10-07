import io
import pytest
from PIL import Image
from nutriomics_final.image_candidates import decode_image


def test_rejects_invalid_and_oversized_uploads():
    for raw in [b'',b'not an image',b'x'*(8*1024*1024+1)]:
        with pytest.raises(ValueError):decode_image(raw)
    buffer=io.BytesIO();Image.new('RGB',(30,20)).save(buffer,format='PNG')
    assert decode_image(buffer.getvalue()).size==(30,20)

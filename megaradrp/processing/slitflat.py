import logging
from numina.flow.processing import TagOptionalCorrector, TagFits

_logger = logging.getLogger('numina.processing')

class SlitFlatCorrector(TagOptionalCorrector):
    '''A Node that corrects a frame from slit flat.'''

    def __init__(self, slitflat, mark=True, tagger=None,datamodel=None,
                 dtype='float32'):

        if tagger is None:
            tagger = TagFits('NUM-BPM', 'Badpixel removed with Numina')

        super(SlitFlatCorrector, self).__init__(datamodel, tagger, mark, dtype)

        self.slitflat = slitflat

    def _run(self, img):
        imgid = self.get_imgid(img)

        _logger.debug('correcting slit flat in %s', imgid)

        img[0].data /= self.slitflat

        return img
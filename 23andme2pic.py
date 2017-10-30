from collections import defaultdict
from gzip import GzipFile
from time import time

import os

from PIL import Image
from PIL.ImageColor import getrgb
from PIL import ImageDraw
from PIL import ImageFilter

from utils import make_search_key

# some constants from VCFUtils
CHROM_PREFIX_VAR_1 = 'chr'
CHROM_PREFIX_VAR_1_LEN = len(CHROM_PREFIX_VAR_1)
CHROM_PREFIX_VAR_2 = 'NC_'
CHROM_PREFIX_VAR_2_LEN = len(CHROM_PREFIX_VAR_2)
X_CHROM_NAME = 'X'
X_CHROM_NUMBER = 23
Y_CHROM_NAME = 'Y'
Y_CHROM_NUMBER = 24


curr_dir = os.path.dirname(__file__)

FREQ_BASE_NAME = "actual_exac_frq_balanced.txt.gz"
REF_BASE_NAME = "ref_base_sorted.txt.gz"

REF_BASES_DIR = "ref_bases"
BACKGROUNDS_DIR = "backgrounds"

EXAC_FREQ_BASE_PATH = os.path.join(curr_dir, REF_BASES_DIR, FREQ_BASE_NAME)
REF_BASE_PATH = os.path.join(curr_dir, REF_BASES_DIR, REF_BASE_NAME)
SAMPLES_FOLDER = os.path.join(curr_dir, "samples")


SOME_SAMPLES = (
    "genome_mario_ricci_v3_Full_20160712125931.txt.gz",
    "genome_none_none_v4_Full_20170626140415.txt.gz",
    "genome_v5_Full_20170826191000.txt.gz",
)

MINIMAL_WIDTH = 600
MINIMAL_HEIGHT = 400
DEFAULT_BORDER_X = 24
DEFAULT_BORDER_Y = 24
DEFAULT_THICKNESS = 2
DEFAULT_IMAGE_MODE = "RGBA"
DEFAULT_EXTENSION = ".png"


WHITE = getrgb("white")
BLACK = getrgb("black")
BLUE = getrgb("blue")
LIGHTBLUE = getrgb("lightblue")
SKYBLUE = getrgb("skyblue")
GREEN = getrgb("green")
GRAY = getrgb("gray")

RED = getrgb("red")
TOMATO = getrgb("tomato")
INDIANRED = getrgb("indianred")
MEDIUMVIOLETRED = getrgb("mediumvioletred")
ORANGE = getrgb("orange")
ORANGERED = getrgb("orangered")
DARKRED = getrgb("darkred")

VIOLET = getrgb("violet")
BLUEVIOLET = getrgb("blueviolet")

YELLOW = getrgb("yellow")
GOLD = getrgb("gold")

MY_YELLOW = getrgb("#fdf53a")
MY_RED = getrgb("#fa1b09")
MY_CYAN = getrgb("#3995ba")
MY_BLUE = getrgb("#16348e")

JUST_ONE_POINT = 0
SIMPLE_CIRCLE = 1
GRADIENT_CIRCLE = 2
FOUR_EDGE_SLIM_STAR = 3


class ExacFreqBase(object):

    def __init__(self, freq_data_path):
        self._data = {}
        iostream = open(freq_data_path, 'rb')

        with GzipFile(fileobj=iostream) as f:
            lines = (line.split() for line in f.readlines())
            for (search_key, freq) in lines:
                self._data[int(search_key)] = float(freq)

    def __getitem__(self, item):
        return self._data[item]

    def get(self, search_key, default=0):
        return self._data.get(search_key, default)


class _ChromData(object):

    def __init__(self):
        self._positions = {}
        self._max_position = None

    @property
    def max_position(self):
        if self._max_position is None:
            self._max_position = max(self._positions.keys())
        return self._max_position

    @property
    def length(self):
        return self.max_position - 1

    def __getitem__(self, item):
        return self._positions[item]

    def __setitem__(self, key, value):
        self._max_position = None
        self._positions[key] = value

    def __len__(self):
        return len(self._positions)

    def keys(self):
        return self._positions.keys()


class RefBase(object):

    def __init__(self, ref_data_path):
        self._sorted_chroms = [str(x) for x in xrange(1, 23)]
        self._sorted_chroms += ["X", "Y"]

        self._ref_dict = {chrom: _ChromData() for chrom in self._sorted_chroms}

        iostream = open(ref_data_path, 'rb')
        with GzipFile(fileobj=iostream) as f:
            lines = (line.split() for line in f.readlines())
            for (chrom, pos, _, ref) in lines:
                self._ref_dict[chrom][int(pos)] = ref

        # chroms absolute positions
        self._chroms_abs_pos = {}
        current_pos = 0
        for chrom in self._sorted_chroms:
            self._chroms_abs_pos[chrom] = current_pos
            current_pos += self._ref_dict[chrom].length

        self._size = current_pos

    def chrom_abs_position(self, chrom, pos):
        chrom_position = self._chroms_abs_pos.get(chrom)
        if chrom_position is not None:
            result = chrom_position + pos
        else:
            result = 0

        return result

    def __getitem__(self, item):
        return self._ref_dict[item]

    def all_pos_iter(self):
        index = 0
        for chrom in self.chromosomes:
            sorted_positions = self._ref_dict[chrom].keys()
            sorted_positions.sort()
            for pos in sorted_positions:
                yield (index, chrom, pos)
                index += 1

    @property
    def size(self):
        return self._size

    @property
    def chromosomes(self):
        return self._sorted_chroms


class Mark(object):

    @staticmethod
    def calc_circle_coords(x, y, r):
        return x - r, y - r, x + r, y + r

    @staticmethod
    def _four_edge_star_polygon(x1, y1, x2, y2, thickness):
        mid_x = x1 + (x2 - x1) // 2
        mid_y = y1 + (y2 - y1) // 2
        result = [
            (mid_x, y1),
            (mid_x + thickness, mid_y - thickness),
            (x2, mid_y),
            (mid_x + thickness, mid_y + thickness),
            (mid_x, y2),
            (mid_x - thickness, mid_y + thickness),
            (x1, mid_y),
            (mid_x - thickness, mid_y - thickness),
            (mid_x, y1)
        ]
        return result

    @staticmethod
    def _romb_polygon(x1, y1, x2, y2):
        mid_x = x1 + (x2 - x1) // 2
        mid_y = y1 + (y2 - y1) // 2
        result = [
            (mid_x, y1),
            (x2, mid_y),
            (mid_x, y2),
            (x1, mid_y),
            (mid_x, y1)
        ]
        return result

    @staticmethod
    def _triangle_polygon(x1, y1, x2, y2):
        mid_x = x1 + (x2 - x1) // 2
        result = [
            (mid_x, y1),
            (x2, y2),
            (x1, y2),
            (mid_x, y1),
        ]
        return result

    @staticmethod
    def decrease_color(color, step_r=0, step_g=0, step_b=0, default_step=10):
        color_len = len(color)
        if color_len == 4:
            r, g, b, a = color
        elif color_len == 3:
            r, g, b = color
        else:
            raise ValueError

        step_r = step_r or default_step
        step_g = step_g or default_step
        step_b = step_b or default_step

        r = max(0, r - step_r)
        g = max(0, g - step_g)
        b = max(0, b - step_b)

        if color_len == 4:
            return (r, g, b, a)
        else:
            return (r, g, b)

    @staticmethod
    def increase_color(color, step_r=0, step_g=0, step_b=0, default_step=10):
        color_len = len(color)
        if color_len == 4:
            r, g, b, a = color
        elif color_len == 3:
            r, g, b = color
        else:
            raise ValueError

        step_r = step_r or default_step
        step_g = step_g or default_step
        step_b = step_b or default_step

        r = max(0, r + step_r)
        g = max(0, g + step_g)
        b = max(0, b + step_b)

        if color_len == 4:
            return (r, g, b, a)
        else:
            return (r, g, b)

    def __init__(self, mark_type, r, color, **kwargs):
        self._mark_type = mark_type
        self._r = r
        self._color = color
        self._mode = kwargs.get("mode", DEFAULT_IMAGE_MODE)
        self._thickness = kwargs.get("thickness", DEFAULT_THICKNESS)
        self._outline_color = kwargs.get("outline_color")
        if mark_type == JUST_ONE_POINT:
            self._image = None
        elif mark_type == SIMPLE_CIRCLE:
            self._image = self._make_simple_circle()
        elif mark_type == GRADIENT_CIRCLE:
            self._image = self._make_gradient_circle()
        elif mark_type == FOUR_EDGE_SLIM_STAR:
            self._image = self._make_four_edge_star()
        else:
            raise ValueError("Unknown Mark Type <%d>" % mark_type)

    def _make_simple_circle(self):
        mode = self._mode
        r = self._r
        image = Image.new(
            mode,
            (r * 2 + 1, r * 2 + 1),
        )
        draw = ImageDraw.Draw(image)

        if r == 2:
            polygon = [
                (1, 0),  (3, 0),
                (4, 1),  (4, 3),
                (3, 4),  (1, 4),
                (0, 3),  (0, 1),
                (1, 0),
            ]
            draw.polygon(
                polygon,
                fill=self._color,
                outline=self._outline_color
            )
        else:
            draw.ellipse(
                (0, 0, r * 2 + 1, r * 2 + 1),
                fill=self._color,
                outline=self._outline_color
            )

        return image

    def _make_four_edge_star(self):
        mode = self._mode
        r = self._r
        image = Image.new(
            mode,
            (r * 2 + 1, r * 2 + 1),
        )

        draw = ImageDraw.Draw(image)

        polygon = self._four_edge_star_polygon(
            0,
            0,
            r * 2 + 1,
            r * 2 + 1,
            self._thickness
        )

        draw.polygon(
            polygon,
            fill=self._color,
            outline=self._outline_color
        )

        return image

    def _make_gradient_circle(self):
        start_color = self._color
        mode = self._mode
        r = self._r
        x0, y0 = r + 1, r + 1
        image = Image.new(
            mode,
            (r * 2 + 1, r * 2 + 1),
        )
        draw = ImageDraw.Draw(image)
        start_radius = 2
        draw.ellipse(
            self.calc_circle_coords(x0, y0, start_radius),
            start_color,
            start_color
        )

        radius = start_radius
        color = start_color

        for _ in range(r - start_radius - 1):
            radius += 1
            color = self.decrease_color(color)
            draw.ellipse(
                self.calc_circle_coords(x0, y0, radius),
                outline=color,
            )

        return image

    def _make_gradient_romb(self):
        # Under construction...
        start_color = self._color
        mode = self._mode
        r = self._r
        x0, y0 = r + 1, r + 1
        image = Image.new(
            mode,
            (r * 2 + 1, r * 2 + 1),
        )
        draw = ImageDraw.Draw(image)
        start_radius = 2

        return image

    @property
    def image(self):
        return self._image

    @property
    def color(self):
        return self._color

    @property
    def radius(self):
        return self._r

    @property
    def width(self):
        return self._image.width

    @property
    def height(self):
        return self._image.height

    @property
    def type(self):
        return self._mark_type


class MarkGroup(object):

    """Base mark group class."""

    ALLELES = None
    COLOR = None
    OUTLINE_COLOR = None

    def __init__(self, alleles=None, color=None, **kwargs):
        self._alleles = alleles or self.ALLELES
        if not isinstance(self._alleles, (tuple, list)):
            self._alleles = (self._alleles,)
        color = color or self.COLOR
        self._color = getrgb(color) if isinstance(color, basestring) else color

        outline_color = kwargs.get("outline_color", self.OUTLINE_COLOR)
        if isinstance(outline_color, basestring):
            self._outline_color = getrgb(outline_color)
        else:
            self._outline_color = outline_color

        self._images = {
            1: Mark(
                FOUR_EDGE_SLIM_STAR, 30,
                color=self._color,
                thicknes=3,
                outline_color=self._outline_color
            ).image,
            2: Mark(
                FOUR_EDGE_SLIM_STAR, 18,
                color=self._color,
                thicknes=2,
                outline_color=self._outline_color
            ).image,
            3: Mark(
                FOUR_EDGE_SLIM_STAR, 12,
                color=self._color,
                thicknes=1,
                outline_color=self._outline_color
            ).image,
            4: Mark(
                FOUR_EDGE_SLIM_STAR, 6,
                color=self._color,
                thicknes=1,
                outline_color=self._outline_color
            ).image,
            5: Mark(
                SIMPLE_CIRCLE, 2,
                color=self._color,
                outline_color=self._outline_color
            ).image,
            6: Mark(
                SIMPLE_CIRCLE, 2,
                color=WHITE
            ).image,
            7: Mark(
                SIMPLE_CIRCLE, 1,
                color=WHITE
            ).image,
            8: Mark(
                JUST_ONE_POINT, 1,
                color=WHITE
            )
        }

    def image(self, freq_value):
        return self._images[freq_value]

    @property
    def alleles(self):
        return self._alleles

    @property
    def color(self):
        return self._color


class AA(MarkGroup):

    ALLELES = ("A",)
    COLOR = WHITE
    OUTLINE_COLOR = MY_RED


class CC(MarkGroup):

    ALLELES = ("C",)
    OUTLINE_COLOR = MY_YELLOW
    COLOR = WHITE


class GG(MarkGroup):

    ALLELES = ("G",)
    COLOR = WHITE
    OUTLINE_COLOR = MY_CYAN


class TT(MarkGroup):

    ALLELES = ("T",)
    COLOR = WHITE
    OUTLINE_COLOR = MY_BLUE


class CategoriesMap(object):

    DEFAULT_CATEGORIES_SIZES = (
        10, 20, 50, 100, 200, 500, 750, 750
    )

    MINIMAL_CATEGORIES_COUNT = 5

    def __init__(self, categories_sizes_tuple=None):
        self._categories_sizes = (
            categories_sizes_tuple or self.DEFAULT_CATEGORIES_SIZES
        )
        assert len(self._categories_sizes) > self.MINIMAL_CATEGORIES_COUNT

    def __len__(self):
        return len(self._categories_sizes)

    def categories(self, ratio):
        assert ratio > 0
        result = []
        limit = 0
        for cat_size in self._categories_sizes:
            result.append((limit, int(limit + cat_size * ratio)))
            limit += int(cat_size * ratio)

        print result
        return result


class Imaginator(object):

    mark_groups = (
        AA(),
        CC(),
        GG(),
        TT(),
    )

    DEFAULT_COLOR = BLACK

    def __init__(self, ref_base_path, freq_base_path, **kwargs):
        self._rb = RefBase(ref_base_path)
        self._fb = ExacFreqBase(freq_base_path or EXAC_FREQ_BASE_PATH)
        self._categories_map = CategoriesMap()

    def _weed_sample_data(self, sample_data, ratio=1):
        categories = []
        sample_data_list = []
        for chrom in sample_data.keys():
            for pos, (alt, freq) in sample_data[chrom].items():
                sample_data_list.append((chrom, pos, alt, freq))

        sample_data_list = sorted(sample_data_list, key=lambda t: t[3])

        for cat_num, (range_start, range_stop) in enumerate(
                self._categories_map.categories(ratio), 1
        ):
            categories.append(
                (cat_num, sample_data_list[range_start:range_stop])
            )

        categories.reverse()

        return categories

    @staticmethod
    def _get_alt_and_genotype(ref, alleles):
        a = alleles[0]

        if len(alleles) > 1:
            b = alleles[1]
            if a != b:
                if ref == a:
                    alt = b
                    genotype = "0/1"
                elif ref == b:
                    alt = a
                    genotype = "0/1"
                else:
                    alt = "%s,%s" % (a, b)
                    genotype = "1/2"
            else:
                if ref == a:
                    alt = "."
                    genotype = "0/0"
                else:
                    alt = a
                    genotype = "1/1"
        else:
            if ref == a:
                alt = "."
                genotype = "0"
            else:
                alt = a
                genotype = "1"

        return alt, genotype

    @staticmethod
    def _parse_chrom_column(data):
        if data.startswith(CHROM_PREFIX_VAR_1):
            data = data[CHROM_PREFIX_VAR_1_LEN:]
        elif data.startswith(CHROM_PREFIX_VAR_2):
            data = data[CHROM_PREFIX_VAR_2_LEN:]
            dot_index = data.find('.')
            if dot_index != -1:
                data = data[:dot_index]

            return int(data)

        if data == X_CHROM_NAME:
            return X_CHROM_NUMBER
        if data == Y_CHROM_NAME:
            return Y_CHROM_NUMBER

        return int(data)

    def _read_sample_data(self, sample_path):
        if not os.path.exists(sample_path):
            raise RuntimeError("Sample path %s does not exist" % sample_path)

        sample_data = defaultdict(lambda: defaultdict(dict))

        if sample_path.endswith(".gz"):
            input_stream = open(sample_path, 'rb')
            input_stream = GzipFile(fileobj=input_stream)
        else:
            input_stream = open(sample_path, 'r')

        try:
            for line in input_stream:
                if line.startswith("#"):
                    continue

                rsid, chrom, pos, alleles = line.split()

                if rsid.startswith("i"):
                    continue

                # ignore mitochondrial chromosomes
                if chrom.startswith("M"):
                    continue

                alleles = alleles.rstrip().upper()

                # ignore "indels"
                if any([a not in "ACGT" for a in alleles]):
                    continue

                pos = int(pos)

                try:
                    ref = self._rb[chrom][pos]
                except KeyError:
                    continue

                alt, genotype = self._get_alt_and_genotype(ref, alleles)

                if alt == ".":
                    # store just mutations, no matches
                    continue

                search_key = make_search_key(
                    self._parse_chrom_column(chrom),
                    pos, ref, alt
                )

                freq = self._fb.get(search_key)
                if not freq:
                    continue

                sample_data[chrom][pos] = (alt, freq)

        finally:
            input_stream.close()

        return sample_data

    def _make_image(self, sample_data, **kwargs):
        border_x = kwargs.get("border_x", DEFAULT_BORDER_X)
        border_y = kwargs.get("border_y", DEFAULT_BORDER_Y)
        background_image = kwargs.get("background_image")

        if background_image:
            result_image = Image.open(background_image)
            image_width = result_image.width
            if image_width < MINIMAL_WIDTH:
                raise ValueError(
                    "Image width must be great than %d" % MINIMAL_WIDTH
                )
            image_height = result_image.height
            if image_height < MINIMAL_HEIGHT:
                raise ValueError(
                    "Image height must be great than %d" % MINIMAL_HEIGHT
                )
        else:
            image_width = kwargs.get("width")
            image_height = kwargs.get("height")
            if image_width is None:
                raise ValueError("expected image width")
            if image_height is None:
                raise ValueError("expected image height")
            image_mode = kwargs.get("image_mode", DEFAULT_IMAGE_MODE)
            result_image = Image.new(
                image_mode, (image_width, image_height), BLACK
            )

        internal_filter = kwargs.get("internal_filter")
        internal_filter_repeats = kwargs.get("internal_filter_repeats", 1)
        final_filter = kwargs.get("final_filter")

        chrom_image_width = image_width - 2 * border_x
        chrom_image_height = image_height - 2 * border_y

        chrom_image_size = chrom_image_width * chrom_image_height

        ratio = self._rb.size / float(chrom_image_size)

        sample_data = self._weed_sample_data(sample_data, 1234 / ratio)

        pixels = result_image.load()

        # most 3 frequent alleles categories
        for cat_num, cat_data in sample_data[:3]:
            for chrom, pos, alt, _ in cat_data:
                if not alt:
                    continue

                mark_group = None
                for group in self.mark_groups:
                    if alt in group.alleles:
                        mark_group = group
                        break

                    if mark_group is None:
                        continue

                mark = mark_group.image(cat_num)
                current_main_position = self._rb.chrom_abs_position(
                    chrom, pos
                )
                current_image_position = int(current_main_position / ratio)
                if current_main_position == 0:
                    continue

                y, x = divmod(current_image_position, chrom_image_width)

                if isinstance(mark, Mark) and mark.type == JUST_ONE_POINT:
                    pixels[border_x + x, border_y + y] = mark.color
                else:
                    result_image.paste(
                        mark,
                        (
                            border_x + x - mark.width // 2,
                            border_y + y - mark.height // 2
                        ),
                        mark
                    )

        if internal_filter:
            for _ in range(internal_filter_repeats):
                result_image = result_image.filter(internal_filter)

        # other categories
        # TODO: Do not Repeat Yourself
        for cat_num, cat_data in sample_data[3:]:
            for chrom, pos, alt, _ in cat_data:
                if not alt:
                    continue

                mark_group = None
                for group in self.mark_groups:
                    if alt in group.alleles:
                        mark_group = group
                        break

                    if mark_group is None:
                        continue

                mark = mark_group.image(cat_num)
                current_main_position = self._rb.chrom_abs_position(
                    chrom, pos
                )
                if current_main_position == 0:
                    continue

                current_image_position = int(current_main_position / ratio)

                y, x = divmod(current_image_position, chrom_image_width)

                if isinstance(mark, Mark) and mark.type == JUST_ONE_POINT:
                    pixels[border_x + x, border_y + y] = mark.color
                else:
                    result_image.paste(
                        mark,
                        (
                            border_x + x - mark.width // 2,
                            border_y + y - mark.height // 2
                        ),
                        mark
                    )

        if final_filter:
            result_image = result_image.filter(final_filter)

        return result_image

    def make(self, sample_path, result_path, **kwargs):
        sample_data = self._read_sample_data(sample_path)
        main_image = self._make_image(sample_data, **kwargs)
        main_image.save(result_path)


imaginator = Imaginator(
    REF_BASE_PATH,
    EXAC_FREQ_BASE_PATH
)

for sample in SOME_SAMPLES:
    sample_path = os.path.join(SAMPLES_FOLDER, sample)
    image_path = os.path.join("/home/anton", sample + ".png")
    t0 = time()

    imaginator.make(
        sample_path, image_path,
        # background_image=os.path.join(curr_dir, BACKGROUNDS_DIR, "scaled_bg04.png"),
        background_image=os.path.join(curr_dir, BACKGROUNDS_DIR, "nightsky_900_600.jpg"),
        # background_image=os.path.join(curr_dir, BACKGROUNDS_DIR, "nightsky_4275_2850.jpg"),
        width=800,
        height=600,
        border_x=36,
        border_y=36,
        internal_filter=ImageFilter.SMOOTH_MORE,
        internal_filter_repeats=3,
        final_filter=ImageFilter.SMOOTH
    )
    print "image for <%s> has been created (dt = %s)" % (sample, time() - t0)

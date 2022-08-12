from gfftools.helper import _get, _get_and_cast


def test_get_with_different_iteratbles():

    l = ['a', 'b', 'c']
    t = ('a', 'b', 'c')

    for x in [l,t]:
        assert _get(1, x, None) == 'b'
        assert _get(4, x, None) == None
        assert _get(-1000, x, None) == None


class TestGetAndCast:

    def test_get_and_cast_with_correct_list(self):
        l = ['1', '2', '3']

        assert _get_and_cast(1, l, int, None) == 2
        assert _get_and_cast(7, l, int, None) == None
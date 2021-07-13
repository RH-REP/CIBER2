from . import shelf

class Child(shelf.Reduction_Shelf):
    def test(self):
        print("hello")
    # def bye(self):　#子クラスで新たに定義したメソッド
    #     print("good bye")
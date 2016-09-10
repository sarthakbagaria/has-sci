# Has-Sci

A collection of computational methods in science.

The aim of this libary is to utilise the powerful type system of haskell and its awesome expressivity to keep the underlying mathematics visible, while making extensive use of typeclasses to allow easy plugging of external objects, methods and numerical libraries.


## Simulators

Build and run using [stack](https://docs.haskellstack.org/en/stable/README/#how-to-install):
```
stack build
stack exec has-sci-physics-examples-exe
```

For some simulators the initial simulation window may be tiny depending on the configuration. You can resize the window and use Ctrl + Drag-Down to zoom in or Ctrl + Drag-Up to zoom out.

You can change the simulator from the [main executable source](https://github.com/sarthakbagaria/has-sci/blob/master/has-sci-physics-examples/app/Main.hs).
And you can change configuration of a simulator from its source module.
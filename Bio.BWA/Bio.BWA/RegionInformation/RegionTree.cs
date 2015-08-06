using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace Bio.BWA
{
    /// <summary>
    ///  Holds genomic regions based on Arne Andersson Self Balancing Binary Search Tree.
    /// </summary>
    public class RegionTree
    {
        #region AATreeNode Class
        [DebuggerDisplay("Value: {Value}")]
        internal class AATreeNode
        {
            /// <summary>
            /// Constructor to initialize null node.
            /// </summary>
            public AATreeNode()
            {
                this.Left = this;
                this.Right = this;
                this.Level = 0;
            }

            public AATreeNode(Region value, AATreeNode nullNode)
            {
                this.Value = value;
                this.Left = nullNode;
                this.Right = nullNode;
                this.Level = 1;
            }

            public Region Value { get; set; }
            public AATreeNode Left { get; set; }
            public AATreeNode Right { get; set; }
            public int Level { get; set; }
        }
        #endregion

        #region Member variables
        /// <summary>
        /// Holds null node.
        /// </summary>
        private static AATreeNode NullNode = new AATreeNode(new Region(String.Empty, 0, 1), RegionTree.NullNode);

        /// <summary>
        /// Holds Root node.n
        /// </summary>
        private AATreeNode root;

        #endregion

        #region Constructors

        /// <summary>
        /// Initializes an instance of AATree class with specified comparer.
        /// </summary>
        public RegionTree()
        {
            this.root = NullNode;
            this.DefaultValue = default(Region);
        }
        #endregion

        #region Properties
        /// <summary>
        /// Gets or sets Default value.
        /// By default this is set to default value of T.
        /// </summary>
        public Region DefaultValue { get; set; }

        /// <summary>
        /// Gets number of elements present in the AATree.
        /// </summary>
        public long Count { get; private set; }
        #endregion

        #region Methods
        /// <summary>
        /// Tries to add specified value to the AATree.
        /// If the value is already present in the tree then this method returns without adding.
        /// </summary>
        /// <param name="value">Value to add.</param>
        /// <returns>Returns true if value is added successfully, else returns false.</returns>
        public void Add(Region value)
        {
            bool result = false;
            if (this.root == RegionTree.NullNode)
            {
                this.root = new AATreeNode(value, RegionTree.NullNode);
                result = true;
            }

            if (!result)
            {
                AATreeNode node = this.root;
                Stack<AATreeNode> parentNodes = new Stack<AATreeNode>();
                while (true)
                {
                    int keyCompResult = value.CompareTo(node.Value);
                    if (keyCompResult == 0)
                    {
                        // key already exists.
                        result = false;
                        node.Value.Merge (value);
                        break;
                    }
                    else if (keyCompResult < 0)
                    {
                        // go to left.
                        if (node.Left == RegionTree.NullNode)
                        {
                            node.Left = new AATreeNode(value, RegionTree.NullNode);
                            result = true;
                            break;
                        }
                        else
                        {
                            parentNodes.Push(node);
                            node = node.Left;
                        }
                    }
                    else
                    {
                        // go to right.
                        if (node.Right == RegionTree.NullNode)
                        {
                            node.Right = new AATreeNode(value, RegionTree.NullNode);
                            result = true;
                            break;
                        }
                        else
                        {
                            parentNodes.Push(node);
                            node = node.Right;
                        }
                    }
                }

                if (result)
                {
                    parentNodes.Push(node);
                    BalanceTreeAfterAdd(parentNodes);
                }
            }

            if (result)
            {
                Count++;
            }
        }


        /// <summary>
        /// Searches for the specified value in the AATree.
        /// If found returns the node containing the value in node out param, else this param contains NullNode.
        /// </summary>
        /// <param name="value">Value to search.</param>
        /// <param name="node">AATree node.</param>
        /// <returns>Returns true if the value is found, else returns false.</returns>
        internal bool TrySearch(Region value, out Region output)
        {
            bool result = false;
            AATreeNode node;
            node = RegionTree.NullNode;

            var currentNode = this.root;
            while (currentNode != RegionTree.NullNode)
            {
                int compResult = value.CompareTo(currentNode.Value);
                if (compResult == 0)
                {
                    node = currentNode;
                    result = true;
                    break;
                }

                if (compResult < 0)
                {
                    //goto left
                    currentNode = currentNode.Left;
                }
                else
                {
                    //goto right
                    currentNode = currentNode.Right;
                }
            }
            output = node.Value;
            return result;
        }

        /// <summary>
        /// Gets values using InOrder traversal.
        /// </summary>
        public IEnumerable<Region> InOrderTraversal()
        {
            return RegionTree.InOrderTraversal(this.root);
        }

        /// <summary>
        /// Gets values using PreOrder traversal.
        /// </summary>
        public IEnumerable<Region> PreOrderTraversal()
        {
            return RegionTree.PreOrderTraversal(this.root);
        }

        /// <summary>
        /// Gets values using PostOrder traversal.
        /// </summary>
        public IEnumerable<Region> PostOrderTraversal()
        {
            return RegionTree.PostOrderTraversal(this.root);
        }


        // InOrder Traversal implementation.
        private static IEnumerable<Region> InOrderTraversal(AATreeNode node)
        {
            if (node == RegionTree.NullNode)
            {
                yield break;
            }

            Stack<AATreeNode> nodes = new Stack<AATreeNode>();
            AATreeNode currentNode = node;
            while (nodes.Count > 0 || currentNode != RegionTree.NullNode)
            {
                if (currentNode != RegionTree.NullNode)
                {
                    nodes.Push(currentNode);
                    currentNode = currentNode.Left;
                }
                else
                {
                    currentNode = nodes.Pop();
                    yield return currentNode.Value;
                    currentNode = currentNode.Right;
                }
            }
        }

        // PreOrder Traversal implementation.
        private static IEnumerable<Region> PreOrderTraversal(AATreeNode node)
        {
            if (node == RegionTree.NullNode)
            {
                yield break;
            }

            Stack<AATreeNode> stack = new Stack<AATreeNode>();
            stack.Push(node);
            AATreeNode currentNode = RegionTree.NullNode;
            while (stack.Count > 0)
            {
                currentNode = stack.Pop();
                if (currentNode.Right != RegionTree.NullNode)
                {
                    stack.Push(currentNode.Right);
                }

                if (currentNode.Left != RegionTree.NullNode)
                {
                    stack.Push(currentNode.Left);
                }

                yield return currentNode.Value;
            }
        }

        // PostOrder Traversal implementation.
        private static IEnumerable<Region> PostOrderTraversal(AATreeNode node)
        {
            if (node == RegionTree.NullNode)
            {
                yield break;
            }

            Stack<AATreeNode> nodes = new Stack<AATreeNode>();
            AATreeNode currentNode = node;
            while (true)
            {
                if (currentNode != RegionTree.NullNode)
                {
                    if (currentNode.Right != RegionTree.NullNode)
                    {
                        nodes.Push(currentNode.Right);
                    }

                    nodes.Push(currentNode);
                    currentNode = currentNode.Left;
                }
                else
                {
                    if (nodes.Count == 0)
                    {
                        break;
                    }

                    currentNode = nodes.Pop();
                    if (currentNode.Right != RegionTree.NullNode && nodes.Count > 0 && nodes.Peek() == currentNode.Right)
                    {
                        nodes.Pop();  // remove right;
                        nodes.Push(currentNode); // push current
                        currentNode = currentNode.Right;
                    }
                    else
                    {
                        yield return currentNode.Value;
                        currentNode = RegionTree.NullNode;
                    }
                }
            }
        }

        // Tree Balancing After Adding an item
        private void BalanceTreeAfterAdd(Stack<AATreeNode> parentNodes)
        {
            AATreeNode parentNode = null;
            while (parentNodes.Count > 0)
            {
                var node = parentNodes.Pop();
                if (parentNodes.Count > 0)
                {
                    parentNode = parentNodes.Peek();
                }
                else
                {
                    parentNode = null;
                }

                // maintain two properties
                // 1. if a node and its left child have same level then rotate right.
                // 2. if a node, its right child and right node's right child have same level then rotate left.
                if (parentNode != null)
                {
                    bool isLeftNode = parentNode.Left == node;
                    RotateRight(parentNode, isLeftNode ? parentNode.Left : parentNode.Right);
                    RotateLeft(parentNode, isLeftNode ? parentNode.Left : parentNode.Right);
                }
                else
                {
                    // check root node.
                    RotateRight(null, this.root);
                    RotateLeft(null, this.root);
                }
            }
        }

        // Tree Balancing After remving an item
        private void BalanceTreeAfterRemove(Stack<AATreeNode> parentNodes)
        {
            AATreeNode parentNode = null;
            // balance the tree.
            while (parentNodes.Count > 0)
            {
                var node = parentNodes.Pop();
                if ((node.Level - node.Left.Level) > 1 || (node.Level - node.Right.Level) > 1)
                {
                    node.Level--;
                    if (node.Right.Level > node.Level)
                    {
                        node.Right.Level = node.Level;
                    }

                    if (parentNodes.Count > 0)
                    {
                        parentNode = parentNodes.Peek();
                    }
                    else
                    {
                        parentNode = null;
                    }

                    if (parentNode != null)
                    {
                        bool isLeftNode = parentNode.Left == node;

                        if (RotateRight(parentNode, node))
                        {
                            // get the new child from parent.
                            node = isLeftNode ? parentNode.Left : parentNode.Right;
                        }

                        if (RotateRight(node, node.Right))

                        if (RotateRight(node.Right, node.Right.Right))

                        if (RotateLeft(parentNode, node))
                        {
                            node = isLeftNode ? parentNode.Left : parentNode.Right;
                        }

                        RotateLeft(node, node.Right);
                    }
                    else
                    {
                        RotateRight(null, this.root);
                        RotateRight(this.root, this.root.Right);
                        RotateRight(this.root.Right, this.root.Right.Right);
                        RotateLeft(null, this.root);
                        RotateLeft(this.root, this.root.Right);
                    }
                }
            }
        }

        /// <summary>
        /// Split or Rotate left.
        /// </summary>
        /// <param name="parentNode"></param>
        /// <param name="node"></param>
        /// <returns></returns>
        private bool RotateLeft(AATreeNode parentNode, AATreeNode node)
        {
            bool result = false;
            if (node == RegionTree.NullNode)
            {
                return result;
            }

            if (node.Level == node.Right.Level && node.Right.Level == node.Right.Right.Level)
            {
                // rotate left.
                var nodeToMoveUp = node.Right;
                node.Right = nodeToMoveUp.Left;
                nodeToMoveUp.Left = node;
                if (parentNode != null)
                {
                    if (parentNode.Left == node)
                    {
                        parentNode.Left = nodeToMoveUp;
                        parentNode.Left.Level++;
                    }
                    else
                    {
                        parentNode.Right = nodeToMoveUp;
                        parentNode.Right.Level++;
                    }
                }
                else
                {
                    this.root = nodeToMoveUp;
                    this.root.Level++;
                }

                result = true;
            }

            return result;
        }

        /// <summary>
        /// Skew or Rotate right.
        /// </summary>
        /// <param name="parentNode"></param>
        /// <param name="node"></param>
        /// <returns></returns>
        private bool RotateRight(AATreeNode parentNode, AATreeNode node)
        {
            bool result = false;
            if (node == RegionTree.NullNode)
            {
                return result;
            }

            if (node.Level == node.Left.Level)
            {
                // rotate right.
                var nodeToMoveUp = node.Left;
                node.Left = nodeToMoveUp.Right;
                nodeToMoveUp.Right = node;
                if (parentNode != null)
                {
                    if (parentNode.Left == node)
                    {
                        parentNode.Left = nodeToMoveUp;

                    }
                    else
                    {
                        parentNode.Right = nodeToMoveUp;
                    }
                }
                else
                {
                    this.root = nodeToMoveUp;
                }

                result = true;
            }

            return result;
        }

        #endregion
    }


}